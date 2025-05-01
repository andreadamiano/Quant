#ifndef CALL_OPTION_H
#define CALL_OPTION_H

#include "Priceable.h"
#include <memory>

class CallOption : public Priceable {
private:
    double strike;
    double maturity;
    CallOption(double strike, double maturity);  // private constructor

public:
    double price(const BlackScholesModel& model) override;
    static std::unique_ptr<CallOption> newInstance(const double strike, const double maturity);
};

#endif