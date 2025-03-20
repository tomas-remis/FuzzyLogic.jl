#import all required libraies
import Pkg;
Pkg.add("FuzzyLogic");
Pkg.add("Dictionaries");

using FuzzyLogic
using Dictionaries
import FuzzyLogic: FuzzyOr, FuzzyAnd, FuzzyRelation, FuzzyNegation, FuzzyRule

using Plots
#=
Fuzzylogic.jl standard format, input parameter for the fuzzylogic system,
Method of calling the function: fis(service=2, food=3)
=#

fis = @sugfis function traffic_light(queue, density, waiting_T)::green_T
    # Vehicle queue length (unit: meters)
    queue := begin
        domain = 0:100
        short = GaussianMF(0.0, 25.0)      # Short queue (0-50 meters)
        medium = GaussianMF(50.0, 30.0)    # Medium queue (20-80 meters)
        long = GaussianMF(100.0, 35.0)     # Long queue (60-100 meters)
    end

    # Pedestrian density (unit: people per square meter)
    density := begin
        domain = 0:10
        low = TrapezoidalMF(-1, 0, 2, 4)   # Low density (0-4 people/m²)
        medium = TrapezoidalMF(3, 4, 6, 7) # Medium density (4-7 people/m²)
        high = TrapezoidalMF(6, 8, 10, 11) # High density (8-10 people/m²)
    end

    # Pedestrian waiting time (unit: seconds)
    waiting_T := begin
        domain = 0:60
        short = TrapezoidalMF(-5, 0, 10, 20)  # Short wait (0-20 seconds)
        medium = TrapezoidalMF(15, 25, 35, 45) # Medium wait (20-45 seconds)
        long = TrapezoidalMF(40, 50, 60, 65)   # Long wait (50-60 seconds)
    end

    # Vehicle green light duration (unit: seconds)
    green_T := begin
        domain = 0:60
        short = 10
        normal = 30
        long = 2queue, -0.5density, 1.5waiting_T, 15.0
    end

    # 4 optimized fuzzy rules
    queue == long || density == low --> green_T == long
    queue == short && density == high --> green_T == short
    queue == medium && density == medium && waiting_T == medium --> green_T == normal
    waiting_T == long || queue == short --> green_T == short
end

# get the handle of Plots.plot
p = plot(fis)

# save p as PNG image format
savefig(p, "plot.png")

# Membership function expression for C code, need to add more
function to_c(mf::GaussianMF)
    """
    double GaussianMF(double x, double mean, double sigma) {
        return exp(-0.5 * pow((x - mean) / sigma, 2));
    }
    """
end

function to_c(mf::TrapezoidalMF)
    """
    double TrapezoidalMF(double x, double a, double b, double c, double d) {
        if (x <= a || x >= d) return 0.0;
        if (x >= b && x <= c) return 1.0;
        if (x > a && x < b) return (x - a) / (b - a);
        if (x > c && x < d) return (d - x) / (d - c);
        return 0.0;
        }
    """
end

#Generate membership function which is used in code
function generate_mf_definitions(fis::SugenoFuzzySystem)
    visited = DataType[]
    res = ""
    for (var_name, var) in pairs(fis.inputs)
        for (mf_name, mf) in pairs(var.mfs)
            if !(typeof(mf) in visited)
                res *= to_c(mf) * "\n\n"
                push!(visited, typeof(mf))
            end
        end
    end
    return res
end

#Generate Fuzzification part
function collect_properties(x)
    join([getproperty(x, p) for p in propertynames(x)], ", ")
end

function generate_fuzzification(fis::SugenoFuzzySystem)
    res = ""
    for (var_name, var) in pairs(fis.inputs)
        for (mf_name, mf) in pairs(var.mfs)
            line = "\tdouble $(var_name)_$mf_name = $(nameof(typeof(mf)))($var_name, $(collect_properties(mf)));"
            res *= line * "\n"
        end
    end
    return res
end

# generate rule evaluation part 1
function generate_rule_expression(r::FuzzyRelation)
    prop = r.prop
    subj = r.subj
    return "$(subj)_$prop"
end

function generate_rule_expression(r::FuzzyAnd)
    left = generate_rule_expression(r.left)
    right = generate_rule_expression(r.right)
    return "fmin($left, $right)"
end

function generate_rule_expression(r::FuzzyOr)
    left = generate_rule_expression(r.left)
    right = generate_rule_expression(r.right)
    return "fmax($left, $right)"
end

function generate_rule_expression(r::FuzzyNegation)
    prop = r.prop
    subj = r.subj
    return "1-$(subj)_$prop"
end

#Generate rule evaluation part 2
function generate_rules(fis::SugenoFuzzySystem)
    rules_vector = fis.rules
    rules_c_expression = "\n"
    for i in eachindex(rules_vector)
        rule = rules_vector[i]
        ant_return = generate_rule_expression(rule.antecedent)
        rules_c_expression *= "\tdouble rule$i = $ant_return;\n"
    end
    return rules_c_expression * "\n"
end

function generate_rules_consequent(fis::SugenoFuzzySystem)
    rules_vector = fis.rules
    rules_c_expression = "\n"
    for i in eachindex(rules_vector)
        rule = rules_vector[i]
        cons_return = generate_rule_expression(rule.consequent[1])
        rules_c_expression *= "\tdouble r$(i)_out = $cons_return;\n"
    end
    return rules_c_expression
end

#Generate rule outputs
function generate_rule_expression(r::ConstantSugenoOutput)
    return "$(r.c)"
end

function generate_rule_expression(r::LinearSugenoOutput)
    rule_out = ""
    # combine coefficients and unknowns together to append
    rule_out *= join(["$coeff * $var" for (var, coeff) in pairs(r.coeffs)], " + ") *
                " + $(r.offset)"
    return rule_out
end

function generate_outputs(fis::SugenoFuzzySystem)
    rule_out = ""
    i = 1
    for rule in fis.rules
        var_name = rule.consequent[1].prop
        for (var_na, var) in pairs(fis.outputs)
            rule_return = generate_rule_expression(var.mfs[var_name])
            rule_out *= "\tdouble r$(i)_out = $rule_return;\n"
            i = i + 1
        end
    end
    return rule_out * "\n"
end

#Generate Weighted average calculation
function generate_calculation(fis::SugenoFuzzySystem)
    calculation = ""
    calculation *= "\tdouble numerator = "
    calculation *= join(
        ["(rule$rule_idx * r$(rule_idx)_out)"
         for (rule_idx, rule) in enumerate(fis.rules)],
        " + ")
    calculation *= ";\n"
    calculation *= "\tdouble denominator = "
    calculation *= join(
        ["rule$rule_idx" for (rule_idx, rule) in enumerate(fis.rules)], " + ")
    calculation *= ";\n\n"
    calculation *= "\treturn numerator / denominator;\n"
    return calculation * "\n"
end

# main C_generator function
function generate_tip(fis::SugenoFuzzySystem)
    res = ""
    func_def = generate_mf_definitions(fis)
    top_func_param = join(
        [string("double ", ele) for ele in collect(keys(fis.inputs))], ", ")
    top_func_start = "double $(fis.name)($(top_func_param)){\n"
    top_func_body = generate_fuzzification(fis)
    rules = generate_rules(fis)
    rule_outputs = generate_outputs(fis)
    calculation = generate_calculation(fis)
    top_func_end = "}"
    res = func_def * top_func_start * top_func_body * rules * rule_outputs * calculation *
          top_func_end
    return res
end

c_code = generate_tip(fis)
print(c_code)
write("c_code.c", c_code)
