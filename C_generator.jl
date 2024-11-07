import Pkg;
Pkg.add("FuzzyLogic");
Pkg.add("Dictionaries");

using FuzzyLogic
using Dictionaries
import FuzzyLogic: FuzzyOr, FuzzyAnd, FuzzyRelation, FuzzyNegation, FuzzyRule

fis = @sugfis function tipper(service, food)::tip
    service := begin
        domain = 0:10
        poor = GaussianMF(0.0, 1.5)
        good = GaussianMF(5.0, 1.5)
        excellent = GaussianMF(10.0, 1.5)
    end

    food := begin
        domain = 0:10
        rancid = TrapezoidalMF(-2, 0, 1, 3)
        delicious = TrapezoidalMF(7, 9, 10, 12)
    end

    tip := begin
        domain = 0:30
        cheap = 5.002
        average = 15
        generous = 2service, 0.5food, 5.0
    end

    service == poor || food == rancid --> tip == cheap
    service == good --> tip == average
    service == excellent || food == delicious --> tip == generous
end

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

# generate sugeno outputs C code part
function generate_outputs(fis::SugenoFuzzySystem)
    rule_out = ""
    for (var_name, var) in pairs(fis.outputs)
        i = 1
        for (mf_name, mf) in pairs(var.mfs)
            rule_return = generate_rule_expression(mf)
            rule_out *= "\tdouble r$(i)_out = $rule_return;\n"
            i = i + 1
        end
    end
    return rule_out * "\n"
end

function generate_calculation(fis::SugenoFuzzySystem)
    calculation = ""
    for (var_name, var) in pairs(fis.outputs)
        calculation *= "\tdouble numerator = "
        calculation *= join(
            ["(rule$idx * r$(idx)_out)"
             for (idx, (mf_name, mf)) in enumerate(pairs(var.mfs))],
            " + ")
        calculation *= ";\n"
        calculation *= "\tdouble denominator = "
        calculation *= join(
            ["rule$idx" for (idx, (mf_name, mf)) in enumerate(pairs(var.mfs))], " + ")
        calculation *= ";\n\n"
        calculation *= "\treturn numerator / denominator;\n"
    end
    return calculation * "\n"
end

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

function collect_properties(x)
    join([getproperty(x, p) for p in propertynames(x)], ", ")
end

function generate_fuzzification(fis::SugenoFuzzySystem)
    res = ""
    for (var_name, var) in pairs(fis.inputs)
        for (mf_name, mf) in pairs(var.mfs)
            # @show var_name, mf_name, mf
            # println(mf)
            # println(fieldnames(typeof(mf)))
            line = "\tdouble $(var_name)_$mf_name = $(nameof(typeof(mf)))($var_name, $(collect_properties(mf)));"
            # line = "double food_rancid = GaussianMF(food, 1.0, 2.0);" # TODO: fix me
            res *= line * "\n"
        end
    end
    return res
end

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

# main generator function
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

print(generate_tip(fis))
