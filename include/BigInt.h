//
// Created by Ext Culith on 2025/7/27.
//

#ifndef BIGINT_H
#define BIGINT_H
#include <string>
#include <vector>

class BigInt
{
private:
    typedef int eleType;
    typedef std::vector<eleType> vecType;
    bool isMinus;
    bool thenMinus = false;
    vecType vec;
    class bad_input_string { const char* what() { return "bad string"; } };
    BigInt& minusNoSign(const BigInt& num);
    BigInt& addNoSign(const BigInt& num);
    void syncSize();
    void syncMinus();
    void inverse();
    BigInt& divideNoSign(const BigInt& num);
    bool lessNoSign(const BigInt& num) const;
    bool notMoreNoSign(const BigInt& num) const;
    BigInt& multipleNoSign(const BigInt& num);
    BigInt& modNoSign(const BigInt& num);
    vecType getInverse() const;
public:
    BigInt() : isMinus(false), vec() {  }
    explicit BigInt(std::string str);
    BigInt(const char* str) { this->BigInt::BigInt(std::string(str)); }
    template <typename TN> BigInt(TN num) : isMinus(num < 0)
    {
        if (num < 0) num *= -1;
        while (num != 0)
        {
            vec.push_back(num % 10);
            num /= 10;
        }
    }
    void setStr(std::string str);
    void sync();
    BigInt& addNoFlat(const BigInt& num);
    BigInt& multipleInt(int num);
    BigInt& minusNoFlat(const BigInt& num);
    BigInt operator-() { BigInt tmp = *this; tmp.isMinus = true; return tmp; }
    bool operator<(const BigInt& num);
    BigInt& operator+=(const BigInt& num);
    BigInt& operator-=(const BigInt& num);
    BigInt& operator*=(const BigInt& num);
    BigInt& operator^=(const int& num)
    {
        auto tmp = *this;
        for (int i = 1; i < num; ++i)
            multipleNoSign(tmp);
        if (num % 2)
            isMinus = false;
        return *this;
    }
    BigInt& operator/=(const BigInt& num);
    BigInt& operator%=(const BigInt& num)
    {
        return this->modNoSign(num);
    }
    BigInt& operator+(const BigInt& num)
    {
        BigInt tmp = *this;
        tmp += num;
        return tmp;
    }
    BigInt& operator-(const BigInt& num)
    {
        BigInt tmp = *this;
        tmp -= num;
        return tmp;
    }
    BigInt& operator/(const BigInt& num)
    {
        BigInt tmp = *this;
        tmp /= num;
        return tmp;
    }
    BigInt& operator*(const BigInt& num)
    {
        BigInt tmp = *this;
        tmp *= num;
        return tmp;
    }
    BigInt& operator^(const int num)
    {
        BigInt tmp = *this;
        tmp ^= num;
        return tmp;
    }
    eleType& operator[](const unsigned int& n) { return vec[n]; }
    const vecType& ori_vec() { return vec; }
    std::string ostr() const;
    void reserve(int n) { vec.reserve(n); }
};

#endif //BIGINT_H
