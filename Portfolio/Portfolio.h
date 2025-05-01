#ifndef PORTFOLIO_H
#define PORTFOLIO_H

#include "Priceable.h"
#include "Black&scholes.h"
#include <memory>


class Portfolio: public Priceable
{
    public:
        virtual ~Portfolio() {}

        //interface methods 
        virtual size() =0; //n of items in the portfolio
        virtual int add (const double quantity, std::shared_ptr<Priceable> security) =0; //add a new security to the portfolio and return the index 
        virtual void setQuantity(const int index, const double quantity) =0; //modify qunatity at a given index 
        virtual double price(const BlackScholesModel& model) =0; //price of the portfolio

        //factory method 
        static std::shared_ptr<Portfolio> newinstance(); //create a portoflio instance (wihtout exposing the constructor)
}; 


#endif