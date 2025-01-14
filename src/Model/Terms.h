/**
   The Supporting Hyperplane Optimization Toolkit (SHOT).

   @author Andreas Lundell, Åbo Akademi University

   @section LICENSE
   This software is licensed under the Eclipse Public License 2.0.
   Please see the README and LICENSE files for more information.
*/

#pragma once
#include "../Environment.h"
#include "../Enums.h"
#include "../Structs.h"
#include "../Utilities.h"

#include "Variables.h"

#include "ffunc.hpp"

#include <vector>

namespace SHOT
{

using Interval = mc::Interval;
using IntervalVector = std::vector<Interval>;

class Term
{
public:
    double coefficient;

    std::weak_ptr<Problem> ownerProblem;

    virtual double calculate(const VectorDouble& point) const = 0;

    virtual Interval calculate(const IntervalVector& intervalVector) const = 0;

    virtual Interval getBounds();

    void inline takeOwnership(ProblemPtr owner) { ownerProblem = owner; }

    virtual E_Convexity getConvexity() const = 0;

    virtual E_Monotonicity getMonotonicity() const = 0;
};

class LinearTerm : public Term
{
private:
public:
    VariablePtr variable;

    LinearTerm() = default;

    LinearTerm(double coeff, VariablePtr var)
    {
        coefficient = coeff;
        variable = var;
    }

    inline double calculate(const VectorDouble& point) const override
    {
        double value = coefficient * variable->calculate(point);
        return value;
    }

    inline Interval calculate(const IntervalVector& intervalVector) const override
    {
        Interval value = coefficient * variable->calculate(intervalVector);
        return value;
    }

    E_Convexity getConvexity() const override { return E_Convexity::Linear; };

    E_Monotonicity getMonotonicity() const override
    {
        if(coefficient > 0)
            return (E_Monotonicity::Nondecreasing);
        else if(coefficient < 0)
            return (E_Monotonicity::Nonincreasing);
        else
            return (E_Monotonicity::Constant);
    };
};

using LinearTermPtr = std::shared_ptr<LinearTerm>;

inline std::ostream& operator<<(std::ostream& stream, LinearTermPtr term)
{
    if(term->coefficient == 1.0)
    {
        stream << " +";
    }
    else if(term->coefficient == -1.0)
    {
        stream << " -";
    }
    else if(term->coefficient == 0.0)
    {
        stream << " +0.0*";
    }
    else if(term->coefficient > 0)
    {
        stream << " +" << term->coefficient << '*';
    }
    else
    {
        stream << " " << term->coefficient << '*';
    }

    stream << term->variable->name;
    return stream;
}

template <class T> class Terms : public std::vector<T>
{
protected:
    E_Convexity convexity = E_Convexity::NotSet;
    E_Monotonicity monotonicity = E_Monotonicity::NotSet;

    std::weak_ptr<Problem> ownerProblem;

    virtual void updateConvexity() = 0;

    void updateMonotonicity()
    {
        bool areAllNonincreasing = true;
        bool areAllNondecreasing = true;

        for(auto& TERM : *this)
        {
            auto monotonicity = TERM->getMonotonicity();
            areAllNonincreasing = areAllNonincreasing
                && (monotonicity == E_Monotonicity::Nonincreasing || monotonicity == E_Monotonicity::Constant);
            areAllNondecreasing = areAllNondecreasing
                && (monotonicity == E_Monotonicity::Nondecreasing || monotonicity == E_Monotonicity::Constant);
        }

        if(areAllNonincreasing)
            monotonicity = E_Monotonicity::Nonincreasing;
        else if(areAllNondecreasing)
            monotonicity = E_Monotonicity::Nondecreasing;
        else
            monotonicity = E_Monotonicity::Unknown;
    };

public:
    using std::vector<T>::operator[];

    using std::vector<T>::at;
    using std::vector<T>::begin;
    using std::vector<T>::clear;
    using std::vector<T>::end;
    using std::vector<T>::erase;
    using std::vector<T>::push_back;
    using std::vector<T>::reserve;
    using std::vector<T>::resize;
    using std::vector<T>::size;

    Terms() = default;
    Terms(std::initializer_list<T> terms)
    {
        for(auto& TE : terms)
            (*this).push_back(TE);
    };

    double calculate(const VectorDouble& point) const
    {
        double value = 0.0;
        for(auto& TERM : *this)
        {
            value += TERM->calculate(point);
        }

        return value;
    }

    Interval calculate(const IntervalVector& intervalVector) const
    {
        Interval value = Interval(0.0, 0.0);
        for(auto& TERM : *this)
        {
            value += TERM->calculate(intervalVector);
        }

        return value;
    }

    Interval getBounds() const
    {
        Interval bounds(0.0, 0.0);

        for(auto const& TERM : *this)
        {
            bounds += TERM->getBounds();
        }

        return bounds;
    }

    inline void takeOwnership(ProblemPtr owner)
    {
        ownerProblem = owner;

        for(auto& TERM : *this)
        {
            TERM->takeOwnership(owner);
        }
    }

    inline E_Convexity getConvexity()
    {
        if(convexity == E_Convexity::NotSet)
            updateConvexity();

        return (convexity);
    }

    inline E_Monotonicity getMonotonicity()
    {
        if(monotonicity == E_Monotonicity::NotSet)
            updateMonotonicity();

        return (monotonicity);
    }

    inline bool checkAllForConvexityType(E_Convexity convexityType)
    {
        for(auto& TERM : (*this))
        {
            if(TERM->getConvexity() != convexityType)
                return (false);
        }

        return (true);
    }
};

class LinearTerms : public Terms<LinearTermPtr>
{
private:
    void updateConvexity() override { convexity = E_Convexity::Linear; };

public:
    using std::vector<LinearTermPtr>::operator[];

    using std::vector<LinearTermPtr>::at;
    using std::vector<LinearTermPtr>::begin;
    using std::vector<LinearTermPtr>::clear;
    using std::vector<LinearTermPtr>::end;
    using std::vector<LinearTermPtr>::erase;
    using std::vector<LinearTermPtr>::push_back;
    using std::vector<LinearTermPtr>::reserve;
    using std::vector<LinearTermPtr>::resize;
    using std::vector<LinearTermPtr>::size;

    LinearTerms() = default;

    void add(LinearTermPtr term)
    {
        (*this).push_back(term);
        monotonicity = E_Monotonicity::NotSet;
    }

    void add(LinearTerms terms)
    {
        for(auto& TERM : terms)
        {
            (*this).push_back(TERM);
        }

        if(terms.size() > 0)
        {
            monotonicity = E_Monotonicity::NotSet;
        }
    }

    SparseVariableVector calculateGradient([[maybe_unused]] const VectorDouble& point) const
    {
        SparseVariableVector gradient;

        for(auto& T : (*this))
        {
            if(T->coefficient == 0.0)
                continue;

            auto element = gradient.emplace(T->variable, T->coefficient);
            if(!element.second)
            {
                // Element already exists for the variable
                element.first->second += T->coefficient;
            }
        }

        return gradient;
    };
};

class QuadraticTerm : public Term
{
public:
    VariablePtr firstVariable;
    VariablePtr secondVariable;

    bool isBilinear = false;
    bool isSquare = false;
    bool isBinary = false;

    QuadraticTerm() = default;

    QuadraticTerm(double coeff, VariablePtr variable1, VariablePtr variable2)
    {
        coefficient = coeff;
        firstVariable = variable1;
        secondVariable = variable2;

        if(firstVariable != secondVariable)
        {
            isBilinear = true;
        }
        else
        {
            isSquare = true;
        }

        if(firstVariable->properties.type == E_VariableType::Binary
            && secondVariable->properties.type == E_VariableType::Binary)
        {
            isBinary = true;
        }
    };

    inline double calculate(const VectorDouble& point) const override
    {
        double value = coefficient * firstVariable->calculate(point) * secondVariable->calculate(point);
        return value;
    };

    inline Interval calculate(const IntervalVector& intervalVector) const override
    {
        Interval value
            = coefficient * firstVariable->calculate(intervalVector) * secondVariable->calculate(intervalVector);
        return value;
    }

    E_Convexity getConvexity() const override
    {
        if(firstVariable == secondVariable)
        {
            if(coefficient > 0)
            {
                return (E_Convexity::Convex);
            }
            else if(coefficient < 0)
            {
                return (E_Convexity::Concave);
            }
            else
            {
                return (E_Convexity::Linear);
            }
        }

        return (E_Convexity::Nonconvex);
    }

    E_Monotonicity getMonotonicity() const override
    {
        if(coefficient > 0)
        {
            return (E_Monotonicity::Nondecreasing);
        }
        else if(coefficient < 0)
        {
            return (E_Monotonicity::Nonincreasing);
        }
        else
            return (E_Monotonicity::Constant);
    };
};

using QuadraticTermPtr = std::shared_ptr<QuadraticTerm>;

inline std::ostream& operator<<(std::ostream& stream, QuadraticTermPtr term)
{
    if(term->coefficient == 1.0)
    {
        stream << " +";
    }
    else if(term->coefficient == -1.0)
    {
        stream << " -";
    }
    else if(term->coefficient == 0.0)
    {
        stream << " +0.0*";
    }
    else if(term->coefficient > 0)
    {
        stream << " +" << term->coefficient;
    }
    else
    {
        stream << " " << term->coefficient;
    }

    if(term->firstVariable == term->secondVariable)
        stream << term->firstVariable->name << "^2";
    else
        stream << term->firstVariable->name << '*' << term->secondVariable->name;

    return stream;
}

class QuadraticTerms : public Terms<QuadraticTermPtr>
{
private:
    void updateConvexity() override;

public:
    using std::vector<QuadraticTermPtr>::operator[];

    using std::vector<QuadraticTermPtr>::at;
    using std::vector<QuadraticTermPtr>::begin;
    using std::vector<QuadraticTermPtr>::clear;
    using std::vector<QuadraticTermPtr>::end;
    using std::vector<QuadraticTermPtr>::erase;
    using std::vector<QuadraticTermPtr>::push_back;
    using std::vector<QuadraticTermPtr>::reserve;
    using std::vector<QuadraticTermPtr>::resize;
    using std::vector<QuadraticTermPtr>::size;

    QuadraticTerms() = default;

    void add(QuadraticTermPtr term)
    {
        (*this).push_back(term);
        convexity = E_Convexity::NotSet;
        monotonicity = E_Monotonicity::NotSet;
    }

    void add(QuadraticTerms terms)
    {
        for(auto& TERM : terms)
        {
            (*this).push_back(TERM);
        }

        if(terms.size() > 0)
        {
            convexity = E_Convexity::NotSet;
            monotonicity = E_Monotonicity::NotSet;
        }
    }

    SparseVariableVector calculateGradient(const VectorDouble& point) const
    {
        SparseVariableVector gradient;

        for(auto& T : (*this))
        {
            if(T->coefficient == 0.0)
                continue;

            if(T->firstVariable == T->secondVariable) // variable squared
            {
                auto value = 2 * T->coefficient * point[T->firstVariable->index];
                auto element = gradient.emplace(T->firstVariable, value);

                if(!element.second)
                {
                    // Element already exists for the variable
                    element.first->second += value;
                }
            }
            else
            {
                auto value = T->coefficient * point[T->secondVariable->index];
                auto element = gradient.emplace(T->firstVariable, value);

                if(!element.second)
                {
                    // Element already exists for the variable
                    element.first->second += value;
                }

                value = T->coefficient * point[T->firstVariable->index];

                element = gradient.emplace(T->secondVariable, value);

                if(!element.second)
                {
                    // Element already exists for the variable
                    element.first->second += value;
                }
            }
        }

        return gradient;
    };
};

class MonomialTerm : public Term
{
public:
    Variables variables;

    bool isBilinear;
    bool isSquare;
    bool isBinary;

    MonomialTerm()
    {
        isBilinear = false;
        isSquare = false;
        isBinary = false;
    };

    MonomialTerm(double coeff, Variables vars)
    {
        coefficient = coeff;
        variables = vars;

        isBilinear = false;
        isBinary = true;
        isSquare = false;

        for(auto& V : variables)
        {
            if(V->properties.type != E_VariableType::Binary)
            {
                isBinary = false;
                break;
            }
        }
    };

    // Creates a copy of the term, with variables from destinationProblem
    MonomialTerm(const MonomialTerm* term, ProblemPtr destinationProblem);

    inline double calculate(const VectorDouble& point) const override
    {
        double value = coefficient;

        for(auto& V : variables)
        {
            value *= V->calculate(point);
        }

        return value;
    };

    inline Interval calculate(const IntervalVector& intervalVector) const override
    {
        Interval value(coefficient);

        for(auto& V : variables)
        {
            value *= V->calculate(intervalVector);
        }

        return value;
    }

    inline E_Convexity getConvexity() const override { return E_Convexity::Unknown; };

    inline E_Monotonicity getMonotonicity() const override { return E_Monotonicity::Unknown; };
};

using MonomialTermPtr = std::shared_ptr<MonomialTerm>;

inline std::ostream& operator<<(std::ostream& stream, MonomialTermPtr term)
{
    if(term->coefficient == 1.0)
    {
        stream << " +";
    }
    else if(term->coefficient == -1.0)
    {
        stream << " -";
    }
    else if(term->coefficient == 0.0)
    {
        stream << " +0.0";
    }
    else if(term->coefficient > 0)
    {
        stream << " +" << term->coefficient;
    }
    else
    {
        stream << " " << term->coefficient;
    }

    for(auto& V : term->variables)
    {
        stream << '*' << V->name;
    }

    return stream;
}

class MonomialTerms : public Terms<MonomialTermPtr>
{
private:
    void updateConvexity() override
    {
        auto resultConvexity = E_Convexity::Linear;

        for(auto& T : *this)
            resultConvexity = Utilities::combineConvexity(resultConvexity, T->getConvexity());

        convexity = resultConvexity;
    };

public:
    using std::vector<MonomialTermPtr>::operator[];

    using std::vector<MonomialTermPtr>::at;
    using std::vector<MonomialTermPtr>::begin;
    using std::vector<MonomialTermPtr>::clear;
    using std::vector<MonomialTermPtr>::end;
    using std::vector<MonomialTermPtr>::erase;
    using std::vector<MonomialTermPtr>::push_back;
    using std::vector<MonomialTermPtr>::reserve;
    using std::vector<MonomialTermPtr>::resize;
    using std::vector<MonomialTermPtr>::size;

    MonomialTerms() = default;

    void add(MonomialTermPtr term)
    {
        (*this).push_back(term);
        convexity = E_Convexity::NotSet;
        monotonicity = E_Monotonicity::NotSet;
    }

    void add(MonomialTerms terms)
    {
        for(auto& TERM : terms)
        {
            (*this).push_back(TERM);
        }

        if(terms.size() > 0)
        {
            convexity = E_Convexity::NotSet;
            monotonicity = E_Monotonicity::NotSet;
        }
    }

    SparseVariableVector calculateGradient(const VectorDouble& point) const
    {
        SparseVariableVector gradient;

        for(auto& T : (*this))
        {
            if(T->coefficient == 0.0)
                continue;

            for(auto& V1 : T->variables)
            {
                double value = 1.0;

                for(auto& V2 : T->variables)
                {
                    if(V1 == V2)
                        continue;

                    value *= V2->calculate(point);
                }

                gradient.emplace(V1, value);
            }
        };

        return gradient;
    };

    SparseVariableMatrix calculateHessian(const VectorDouble& point) const
    {
        SparseVariableMatrix hessian;

        for(auto& T : (*this))
        {
            if(T->coefficient == 0)
                continue;

            for(auto& V1 : T->variables)
            {
                for(auto& V2 : T->variables)
                {
                    if(V1->index >= V2->index)
                        continue;

                    double value = T->coefficient;

                    for(auto& V3 : T->variables)
                    {
                        if(V3 == V1 || V3 == V2)
                            continue;

                        value *= V3->calculate(point);
                    }

                    std::pair<VariablePtr, VariablePtr> variablePair = std::make_pair(V1, V2);

                    auto element = hessian.emplace(variablePair, value);

                    if(!element.second)
                    {
                        // Element already exists for the variable
                        element.first->second += value;
                    }
                }
            }
        }

        return hessian;
    };
};

class SignomialElement
{
private:
public:
    VariablePtr variable;
    double power;

    SignomialElement(VariablePtr variable, double power) : variable(variable), power(power){};

    inline double calculate(const VectorDouble& point) const { return pow(variable->calculate(point), power); }

    inline Interval calculate(const IntervalVector& intervalVector) const
    {
        auto variableBound = variable->calculate(intervalVector);

        double intpart;
        bool isInteger = (std::modf(power, &intpart) == 0.0);
        int integerValue = (int)round(intpart);
        bool isEven = (integerValue % 2 == 0);

        if(variableBound.l() <= 0)
        {
            if(!isInteger)
                variableBound.l(SHOT_DBL_EPS);
            else if(isInteger && power < 0)
                variableBound.l(SHOT_DBL_EPS);
        }

        Interval bounds;

        if(isInteger)
            bounds = pow(variableBound, (int)power);
        else
            bounds = pow(variableBound, power);

        if(isInteger && isEven && bounds.l() <= 0.0)
            bounds.l(0.0);

        return (bounds);
    }

    inline Interval getBounds()
    {
        auto variableBound = variable->getBound();

        double intpart;
        bool isInteger = (std::modf(power, &intpart) == 0.0);

        if(isInteger && power > 0 && variableBound.l() < 0.0)
            variableBound.l(0.0);
        else if(!isInteger && variableBound.l() <= 0.0)
            variableBound.l(SHOT_DBL_EPS);
        else if(variableBound.l() <= 0)
            variableBound.l(0.0);

        auto bounds = pow(variableBound, 1.0 / power);

        if(bounds.l() <= 0.0)
            bounds.l(0.0);

        return (bounds);
    }

    inline bool tightenBounds(Interval bound)
    {
        if(bound.l() < 0)
            bound.l(SHOT_DBL_EPS);

        double intpart;
        bool isInteger = (std::modf(power, &intpart) == 0.0);

        if(isInteger && power > 0 && bound.l() < 0.0)
            bound.l(0.0);
        else if(bound.l() <= 0.0)
            bound.l(SHOT_DBL_EPS);
        else if(bound.l() < 0.0)
            bound.l(0.0);

        auto interval = pow(bound, 1.0 / power);

        return (variable->tightenBounds(interval));
    }
};

using SignomialElementPtr = std::shared_ptr<SignomialElement>;
using SignomialElements = std::vector<SignomialElementPtr>;

inline std::ostream& operator<<(std::ostream& stream, SignomialElementPtr element)
{
    if(element->power == 1.0)
        stream << element->variable->name;
    else if(element->power > 0.0)
        stream << element->variable->name << '^' << element->power;
    else
        stream << element->variable->name << "^(" << element->power << ')';

    return stream;
}

class SignomialTerm : public Term
{
public:
    SignomialElements elements;

    SignomialTerm() = default;

    SignomialTerm(double coeff, SignomialElements elems)
    {
        coefficient = coeff;
        elements = elems;
    };

    // Creates a copy of the term, with variables from destinationProblem
    SignomialTerm(const SignomialTerm* term, ProblemPtr destinationProblem);

    inline double calculate(const VectorDouble& point) const override
    {
        double value = coefficient;

        for(auto& E : elements)
        {
            value *= E->calculate(point);
        }

        return value;
    };

    inline Interval calculate(const IntervalVector& intervalVector) const override
    {
        Interval value(coefficient);

        for(auto& E : elements)
        {
            value *= E->calculate(intervalVector);
        }

        return value;
    }

    inline E_Convexity getConvexity() const override
    {
        size_t numberPositivePowers = 0;
        double sumPowers = 0.0;

        for(auto& E : elements)
        {
            if(E->power > 0)
            {
                numberPositivePowers++;
            }

            sumPowers += E->power;
        }

        if(elements.size() == 1 && sumPowers == 1.0)
            return (E_Convexity::Linear);

        if(coefficient > 0)
        {
            if(numberPositivePowers == 1 && sumPowers > 1.0)
                return (E_Convexity::Convex);

            if(elements.size() == 1 && sumPowers > 0.0 && sumPowers < 1.0)
                return (E_Convexity::Concave);

            if(numberPositivePowers == 0)
                return (E_Convexity::Convex);

            return (E_Convexity::Nonconvex);
        }
        else if(coefficient < 0)
        {
            if(numberPositivePowers == 1 && sumPowers > 1.0)
                return (E_Convexity::Concave);

            if(numberPositivePowers == elements.size() && sumPowers > 0.0 && sumPowers <= 1.0)
                return (E_Convexity::Convex);

            if(numberPositivePowers == 0)
                return (E_Convexity::Concave);
        }

        return E_Convexity::Nonconvex;
    };

    inline E_Monotonicity getMonotonicity() const override
    {
        size_t numberPositivePowers = 0;
        double sumPowers = 0.0;

        if(coefficient == 0.0)
            return (E_Monotonicity::Constant);

        for(auto& E : elements)
        {
            if(E->power > 0)
            {
                numberPositivePowers++;
            }

            sumPowers += E->power;
        }

        if(coefficient > 0)
        {
            if(elements.size() == 1 && sumPowers == 0.0)
                return (E_Monotonicity::Constant);

            if(elements.size() == 1 && sumPowers > 0.0)
                return (E_Monotonicity::Nondecreasing);

            if(elements.size() == 1 && sumPowers < 0.0)
                return (E_Monotonicity::Nonincreasing);

            if(numberPositivePowers == 0)
            {
                return (E_Monotonicity::Nonincreasing);
            }

            if(numberPositivePowers == elements.size())
            {
                return (E_Monotonicity::Nondecreasing);
            }
        }
        else if(coefficient < 0)
        {
            if(elements.size() == 1 && sumPowers == 0.0)
                return (E_Monotonicity::Constant);

            if(elements.size() == 1 && sumPowers > 0.0)
                return (E_Monotonicity::Nonincreasing);

            if(elements.size() == 1 && sumPowers < 0.0)
                return (E_Monotonicity::Nondecreasing);

            if(numberPositivePowers == 0)
            {
                return (E_Monotonicity::Nondecreasing);
            }

            if(numberPositivePowers == elements.size())
            {
                return (E_Monotonicity::Nonincreasing);
            }
        }

        return E_Monotonicity::Unknown;
    };
};

using SignomialTermPtr = std::shared_ptr<SignomialTerm>;

inline std::ostream& operator<<(std::ostream& stream, SignomialTermPtr term)
{
    if(term->coefficient == 1.0)
    {
        stream << " +";
    }
    else if(term->coefficient == -1.0)
    {
        stream << " -";
    }
    else if(term->coefficient == 0.0)
    {
        stream << " +0.0";
    }
    else if(term->coefficient > 0)
    {
        stream << " +" << term->coefficient;
    }
    else
    {
        stream << " " << term->coefficient;
    }

    for(auto& E : term->elements)
    {
        stream << '*' << E;
    }

    return stream;
}

class SignomialTerms : public Terms<SignomialTermPtr>
{
private:
    void updateConvexity() override
    {
        auto resultConvexity = E_Convexity::Linear;

        for(auto& T : *this)
            resultConvexity = Utilities::combineConvexity(resultConvexity, T->getConvexity());

        convexity = resultConvexity;
    };

public:
    using std::vector<SignomialTermPtr>::operator[];

    using std::vector<SignomialTermPtr>::at;
    using std::vector<SignomialTermPtr>::begin;
    using std::vector<SignomialTermPtr>::clear;
    using std::vector<SignomialTermPtr>::end;
    using std::vector<SignomialTermPtr>::erase;
    using std::vector<SignomialTermPtr>::push_back;
    using std::vector<SignomialTermPtr>::reserve;
    using std::vector<SignomialTermPtr>::resize;
    using std::vector<SignomialTermPtr>::size;

    SignomialTerms() = default;

    void add(SignomialTermPtr term)
    {
        (*this).push_back(term);
        convexity = E_Convexity::NotSet;
        monotonicity = E_Monotonicity::NotSet;
    }

    void add(SignomialTerms terms)
    {
        for(auto& TERM : terms)
        {
            (*this).push_back(TERM);
        }

        if(terms.size() > 0)
        {
            convexity = E_Convexity::NotSet;
            monotonicity = E_Monotonicity::NotSet;
        }
    }

    inline SparseVariableVector calculateGradient(const VectorDouble& point) const
    {
        SparseVariableVector gradient;

        for(auto& T : (*this))
        {
            if(T->coefficient == 0.0)
                continue;

            for(auto& E1 : T->elements)
            {
                double value = 1.0;

                for(auto& E2 : T->elements)
                {
                    if(E1 == E2)
                    {
                        if(E2->power != 1.0)
                            value *= E2->power * pow(E2->variable->calculate(point), E2->power - 1.0);
                    }
                    else
                    {
                        value *= E2->calculate(point);
                    }
                }

                auto element = gradient.emplace(E1->variable, T->coefficient * value);

                if(!element.second)
                {
                    // Element already exists for the variable
                    element.first->second += value;
                }
            }
        };

        return gradient;
    };

    SparseVariableMatrix calculateHessian(const VectorDouble& point) const
    {
        SparseVariableMatrix hessian;

        for(auto& T : (*this))
        {
            if(T->coefficient == 0)
                continue;

            auto value = T->calculate(point);

            for(auto& E1 : T->elements)
            {
                for(auto& E2 : T->elements)
                {
                    if(E1->variable->index > E2->variable->index)
                        continue;

                    double corrFactor;

                    if(E1->variable->index == E2->variable->index)
                    {
                        corrFactor = E1->power * (E1->power - 1.0)
                            / (E1->variable->calculate(point) * E1->variable->calculate(point));
                    }
                    else
                    {
                        corrFactor
                            = E1->power * E2->power / (E1->variable->calculate(point) * E2->variable->calculate(point));
                    }

                    auto variablePair = std::make_pair(E1->variable, E2->variable);

                    auto element = hessian.emplace(variablePair, corrFactor * value);

                    if(!element.second)
                    {
                        // Element already exists for the variable
                        element.first->second += corrFactor * value;
                    }
                }
            }
        }

        return hessian;
    };
};

inline std::ostream& operator<<(std::ostream& stream, LinearTerms terms)
{
    if(terms.size() == 0)
        return stream;

    stream << ' ' << terms.at(0);

    for(size_t i = 1; i < terms.size(); i++)
    {
        stream << terms.at(i);
    }

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, QuadraticTerms terms)
{
    if(terms.size() == 0)
        return stream;

    stream << ' ' << terms.at(0);

    for(size_t i = 1; i < terms.size(); i++)
    {
        stream << terms.at(i);
    }

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, MonomialTerms terms)
{
    if(terms.size() == 0)
        return stream;

    stream << ' ' << terms.at(0);

    for(size_t i = 1; i < terms.size(); i++)
    {
        stream << terms.at(i);
    }

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, SignomialTerms terms)
{
    if(terms.size() == 0)
        return stream;

    stream << ' ' << terms.at(0);

    for(size_t i = 1; i < terms.size(); i++)
    {
        stream << terms.at(i);
    }

    return stream;
}
} // namespace SHOT