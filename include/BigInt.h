#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <type_traits>
#include <cctype>

namespace betadatastructures
{
    template <typename DigitT = uint32_t, typename DoubleDigitT = uint64_t>
    class BigInt
    {
        static_assert(std::is_integral<DigitT>::value && std::is_unsigned<DigitT>::value, "DigitT must be an unsigned integral type.");
        static_assert(std::is_integral<DoubleDigitT>::value && std::is_unsigned<DoubleDigitT>::value, "DoubleDigitT must be an unsigned integral type.");
        static_assert(sizeof(DoubleDigitT) >= 2 * sizeof(DigitT), "DoubleDigitT must be at least twice the size of DigitT.");

    public:
        using digit_type = DigitT;
        using double_digit_type = DoubleDigitT;

        static constexpr digit_type BASE_DIGITS = 9;
        static constexpr digit_type BASE = 1000000000;

        class LazyProxy;

    private:
        enum class State { Eager, Lazy };

        bool is_negative_ = false;
        std::vector<digit_type> digits_;
        State state_ = State::Eager;

    public:
        BigInt();
        BigInt(long long val);
        explicit BigInt(const std::string& s);

        template <typename T, typename TT> friend bool operator<(const BigInt<T, TT>& a, const BigInt<T, TT>& b);

        template <typename T, typename TT> friend bool operator==(const BigInt<T, TT>& a, const BigInt<T, TT>& b);

        LazyProxy lazy();
        std::string toString() const;
        bool isZero() const;
        explicit operator bool() const;

        BigInt& operator=(const BigInt& other);
        BigInt& operator+=(const BigInt& other);
        BigInt& operator-=(const BigInt& other);
        BigInt& operator*=(const BigInt& other);
        BigInt& operator/=(const BigInt& other);
        BigInt& operator%=(const BigInt& other);

    private:
        void assertEagerState() const;
        void trimLeadingZeros();
        void normalize();

        static int compareMagnitude(const BigInt& a, const BigInt& b);
        static void addMagnitude(BigInt& a, const BigInt& b);
        static void subtractMagnitude(BigInt& a, const BigInt& b);
        static BigInt multiplyMagnitude(const BigInt& a, const BigInt& b);
        static BigInt karatsubaMultiply(const BigInt& a, const BigInt& b);
        static std::pair<BigInt, BigInt> divideMagnitude(const BigInt& a, const BigInt& b);
        static digit_type divideShort(BigInt& a, digit_type d);
        static BigInt multiplyShort(const BigInt& a, digit_type d);

        friend class LazyProxy;

    public:
        class LazyProxy
        {
        private:
            BigInt& ref_;
        public:
            explicit LazyProxy(BigInt& ref) : ref_(ref) { ref_.state_ = State::Lazy; }
            BigInt& finish() { ref_.normalize(); ref_.state_ = State::Eager; return ref_; }
            LazyProxy& add(const BigInt& other)
            {
                other.assertEagerState();
                ref_.digits_.resize(std::max(ref_.digits_.size(), other.digits_.size()), 0);
                for (size_t i = 0; i < other.digits_.size(); ++i) { ref_.digits_[i] += other.digits_[i]; }
                return *this;
            }
            LazyProxy& subtract(const BigInt& other)
            {
                other.assertEagerState();
                ref_.digits_.resize(std::max(ref_.digits_.size(), other.digits_.size()), 0);
                for (size_t i = 0; i < other.digits_.size(); ++i)
                {
                    ref_.digits_[i] -= other.digits_[i];
                }
                return *this;
            }
        };
    };

    template <typename D, typename DD> BigInt<D, DD> operator+(BigInt<D, DD> a, const BigInt<D, DD>& b)
    {
        a += b; return a;
    }
    template <typename D, typename DD> BigInt<D, DD> operator-(BigInt<D, DD> a, const BigInt<D, DD>& b)
    {
        a -= b;
        return a;
    }
    template <typename D, typename DD> BigInt<D, DD> operator*(BigInt<D, DD> a, const BigInt<D, DD>& b)
    {
        a *= b;
        return a;
    }
    template <typename D, typename DD> BigInt<D, DD> operator/(BigInt<D, DD> a, const BigInt<D, DD>& b)
    {
        a /= b;
        return a;
    }
    template <typename D, typename DD> BigInt<D, DD> operator%(BigInt<D, DD> a, const BigInt<D, DD>& b)
    {
        a %= b;
        return a;
    }
    template <typename D, typename DD> bool operator>(const BigInt<D, DD>& a, const BigInt<D, DD>& b)
    {
        return b < a;
    }
    template <typename D, typename DD> bool operator<=(const BigInt<D, DD>& a, const BigInt<D, DD>& b)
    {
        return !(b < a);
    }
    template <typename D, typename DD> bool operator>=(const BigInt<D, DD>& a, const BigInt<D, DD>& b)
    {
        return !(a < b);
    }
    template <typename D, typename DD> bool operator!=(const BigInt<D, DD>& a, const BigInt<D, DD>& b)
    {
        return !(a == b);
    }
    template <typename D, typename DD> std::ostream& operator<<(std::ostream& os, const BigInt<D, DD>& num)
    {
        os << num.toString();
        return os;
    }

     template <typename D, typename DD>
    BigInt<D, DD>::BigInt() : state_(State::Eager)
    {
        digits_.push_back(0);
    }

    template <typename D, typename DD>
    BigInt<D, DD>::BigInt(long long val) : state_(State::Eager)
    {
        if (val == 0)
        {
            digits_.push_back(0);
            return;
        }
        if (val < 0)
        {
            is_negative_ = true;
            val = -val;
        }
        while (val > 0)
        {
            digits_.push_back(val % BASE);
            val /= BASE;
        }
    }

    template <typename D, typename DD>
    BigInt<D, DD>::BigInt(const std::string& s) : state_(State::Eager)
    {
        if (s.empty())
        {
            throw std::invalid_argument("Invalid string for BigInt: empty string");
        }
        int start_pos = (s[0] == '-');
        if (s.length() == start_pos)
        {
            throw std::invalid_argument("Invalid string for BigInt: sign only");
        }
        for (size_t i = start_pos; i < s.length(); ++i)
        {
            if (!std::isdigit(s[i]))
            {
                throw std::invalid_argument("Invalid character in BigInt string: " + s);
            }
        }
        for (int i = s.length(); i > start_pos; i -= BASE_DIGITS)
        {
            int len = std::min(static_cast<int>(BASE_DIGITS), i - start_pos);
            digits_.push_back(std::stoi(s.substr(i - len, len)));
        }
        this->trimLeadingZeros();
        if (start_pos && !this->isZero())
        {
            is_negative_ = true;
        }
    }

    template <typename D, typename DD>
    typename BigInt<D, DD>::LazyProxy BigInt<D, DD>::lazy()
    {
        assertEagerState();
        return LazyProxy(*this);
    }

    template <typename D, typename DD>
    bool BigInt<D, DD>::isZero() const
    {
        return digits_.size() == 1 && digits_[0] == 0;
    }

    template <typename D, typename DD>
    BigInt<D, DD>::operator bool() const
    {
        return !isZero();
    }

    template <typename D, typename DD>
    std::string BigInt<D, DD>::toString() const
    {
        assertEagerState();
        if (isZero())
        {
            return "0";
        }

        std::stringstream ss;
        if (is_negative_)
        {
            ss << "-";
        }

        ss << digits_.back();
        for (int i = digits_.size() - 2; i >= 0; --i)
        {
            ss << std::setw(BASE_DIGITS) << std::setfill('0') << digits_[i];
        }
        return ss.str();
    }

    template <typename D, typename DD>
    BigInt<D, DD>& BigInt<D, DD>::operator=(const BigInt& other)
    {
        if (this == &other)
        {
            return *this;
        }
        assertEagerState();
        other.assertEagerState();
        this->is_negative_ = other.is_negative_;
        this->digits_ = other.digits_;
        return *this;
    }

    template <typename D, typename DD>
    BigInt<D, DD>& BigInt<D, DD>::operator+=(const BigInt& other)
    {
        assertEagerState();
        other.assertEagerState();
        if (is_negative_ == other.is_negative_)
        {
            addMagnitude(*this, other);
        }
        else
        {
            if (compareMagnitude(*this, other) >= 0)
            {
                subtractMagnitude(*this, other);
            }
            else
            {
                BigInt temp = other;
                subtractMagnitude(temp, *this);
                *this = temp;
            }
        }
        return *this;
    }

    template <typename D, typename DD>
    BigInt<D, DD>& BigInt<D, DD>::operator-=(const BigInt& other)
    {
        assertEagerState();
        other.assertEagerState();
        if (is_negative_ != other.is_negative_)
        {
            addMagnitude(*this, other);
        }
        else
        {
            if (compareMagnitude(*this, other) >= 0)
            {
                subtractMagnitude(*this, other);
            }
            else
            {
                BigInt temp = other;
                subtractMagnitude(temp, *this);
                *this = temp;
                this->is_negative_ = !this->is_negative_;
            }
        }
        return *this;
    }

    template <typename D, typename DD>
    BigInt<D, DD>& BigInt<D, DD>::operator*=(const BigInt& other)
    {
        assertEagerState();
        other.assertEagerState();
        bool result_is_negative = !this->isZero() && !other.isZero() && (is_negative_ != other.is_negative_);
        *this = multiplyMagnitude(*this, other);
        is_negative_ = result_is_negative;
        return *this;
    }

    template <typename D, typename DD>
    BigInt<D, DD>& BigInt<D, DD>::operator/=(const BigInt& other)
    {
        assertEagerState();
        other.assertEagerState();
        if (other.isZero())
        {
            throw std::runtime_error("Division by zero");
        }
        bool result_is_negative = !this->isZero() && !other.isZero() && (is_negative_ != other.is_negative_);
        auto div_pair = divideMagnitude(*this, other);
        *this = div_pair.first;
        is_negative_ = result_is_negative;
        return *this;
    }

    template <typename D, typename DD>
    BigInt<D, DD>& BigInt<D, DD>::operator%=(const BigInt& other)
    {
        assertEagerState();
        other.assertEagerState();
        if (other.isZero())
        {
            throw std::runtime_error("Division by zero");
        }
        bool result_is_negative = is_negative_;
        auto div_pair = divideMagnitude(*this, other);
        *this = div_pair.second;
        is_negative_ = result_is_negative;
        return *this;
    }

    template <typename D, typename DD>
    void BigInt<D, DD>::assertEagerState() const
    {
        if (state_ == State::Lazy)
        {
            throw std::logic_error("Operation forbidden on a BigInt in a lazy state. Call .finish() first.");
        }
    }

    template <typename D, typename DD>
    void BigInt<D, DD>::trimLeadingZeros()
    {
        while (digits_.size() > 1 && digits_.back() == 0)
        {
            digits_.pop_back();
        }
        if (isZero())
        {
            is_negative_ = false;
        }
    }

    template <typename D, typename DD>
    void BigInt<D, DD>::normalize()
    {
        using signed_t = std::make_signed_t<double_digit_type>;
        signed_t carry = 0;
        for (size_t i = 0; i < digits_.size() || carry != 0; ++i)
        {
            if (i == digits_.size())
            {
                digits_.push_back(0);
            }
            signed_t current_val = static_cast<signed_t>(digits_[i]) + carry;
            if (current_val < 0)
            {
                signed_t borrow = (-current_val + BASE - 1) / BASE;
                carry = -borrow;
                digits_[i] = static_cast<digit_type>(current_val + borrow * BASE);
            }
            else
            {
                carry = current_val / BASE;
                digits_[i] = static_cast<digit_type>(current_val % BASE);
            }
        }
        this->trimLeadingZeros();
    }

    template <typename D, typename DD>
    int BigInt<D, DD>::compareMagnitude(const BigInt& a, const BigInt& b)
    {
        if (a.digits_.size() != b.digits_.size())
        {
            return a.digits_.size() < b.digits_.size() ? -1 : 1;
        }
        for (int i = a.digits_.size() - 1; i >= 0; --i)
        {
            if (a.digits_[i] != b.digits_[i])
            {
                return a.digits_[i] < b.digits_[i] ? -1 : 1;
            }
        }
        return 0;
    }

    template <typename D, typename DD>
    void BigInt<D, DD>::addMagnitude(BigInt& a, const BigInt& b)
    {
        digit_type carry = 0;
        size_t n = std::max(a.digits_.size(), b.digits_.size());
        a.digits_.resize(n, 0);

        for (size_t i = 0; i < n; ++i)
        {
            double_digit_type sum = a.digits_[i] + carry;
            if (i < b.digits_.size())
            {
                sum += b.digits_[i];
            }
            a.digits_[i] = static_cast<digit_type>(sum % BASE);
            carry = static_cast<digit_type>(sum / BASE);
        }
        if (carry > 0)
        {
            a.digits_.push_back(carry);
        }
    }

    template <typename D, typename DD>
    void BigInt<D, DD>::subtractMagnitude(BigInt& a, const BigInt& b)
    {
        long long borrow = 0;
        size_t n = a.digits_.size();
        for (size_t i = 0; i < n; ++i)
        {
            long long b_digit = (i < b.digits_.size()) ? b.digits_[i] : 0;
            long long diff = static_cast<long long>(a.digits_[i]) - b_digit - borrow;
            if (diff < 0)
            {
                diff += BASE;
                borrow = 1;
            }
            else
            {
                borrow = 0;
            }
            a.digits_[i] = static_cast<digit_type>(diff);
        }
        a.trimLeadingZeros();
    }

    template<typename D, typename DD>
    BigInt<D, DD> BigInt<D, DD>::multiplyMagnitude(const BigInt& a, const BigInt& b)
    {
        return karatsubaMultiply(a, b);
    }

    template<typename D, typename DD>
    BigInt<D, DD> BigInt<D, DD>::karatsubaMultiply(const BigInt& a, const BigInt& b)
    {
        if (a.isZero() || b.isZero())
        {
            return BigInt(0);
        }
        if (a.digits_.size() < 32 || b.digits_.size() < 32)
        {
            BigInt res;
            res.digits_.resize(a.digits_.size() + b.digits_.size(), 0);
            for (size_t i = 0; i < a.digits_.size(); ++i)
            {
                double_digit_type carry = 0;
                for (size_t j = 0; j < b.digits_.size(); ++j)
                {
                    double_digit_type prod = static_cast<double_digit_type>(a.digits_[i]) * b.digits_[j] + res.digits_[i+j] + carry;
                    res.digits_[i+j] = static_cast<digit_type>(prod % BASE);
                    carry = prod / BASE;
                }
                if (carry > 0)
                {
                    res.digits_[i + b.digits_.size()] += carry;
                }
            }
            res.trimLeadingZeros();
            return res;
        }

        size_t m = std::min(a.digits_.size(), b.digits_.size()) / 2;
        BigInt a_low, a_high, b_low, b_high;
        a_low.digits_.assign(a.digits_.begin(), a.digits_.begin() + m);
        a_high.digits_.assign(a.digits_.begin() + m, a.digits_.end());
        b_low.digits_.assign(b.digits_.begin(), b.digits_.begin() + m);
        b_high.digits_.assign(b.digits_.begin() + m, b.digits_.end());

        a_low.trimLeadingZeros(); a_high.trimLeadingZeros();
        b_low.trimLeadingZeros(); b_high.trimLeadingZeros();

        BigInt z0 = karatsubaMultiply(a_low, b_low);
        BigInt z2 = karatsubaMultiply(a_high, b_high);
        BigInt z1 = karatsubaMultiply(a_low + a_high, b_low + b_high);
        z1 -= z0;
        z1 -= z2;

        if (!z2.isZero())
        {
            z2.digits_.insert(z2.digits_.begin(), 2 * m, 0);
        }
        if (!z1.isZero())
        {
            z1.digits_.insert(z1.digits_.begin(), m, 0);
        }

        BigInt result = z0 + z1 + z2;
        result.trimLeadingZeros();
        return result;
    }

    template<typename D, typename DD>
    D BigInt<D, DD>::divideShort(BigInt& a, D d)
    {
        if (d == 0)
        {
            throw std::runtime_error("Division by zero");
        }
        double_digit_type r = 0;
        for (int i = a.digits_.size() - 1; i >= 0; --i)
        {
            r = r * BASE + a.digits_[i];
            a.digits_[i] = static_cast<digit_type>(r / d);
            r %= d;
        }
        a.trimLeadingZeros();
        return static_cast<digit_type>(r);
    }

    template<typename D, typename DD>
    BigInt<D, DD> BigInt<D, DD>::multiplyShort(const BigInt& a, digit_type d)
    {
        if (d == 0 || a.isZero())
        {
            return BigInt(0);
        }
        BigInt res = a;
        res.is_negative_ = false;
        double_digit_type carry = 0;
        for (size_t i = 0; i < res.digits_.size() || carry > 0; ++i)
        {
            if (i == res.digits_.size())
            {
                res.digits_.push_back(0);
            }
            double_digit_type prod = static_cast<double_digit_type>(res.digits_[i]) * d + carry;
            res.digits_[i] = static_cast<digit_type>(prod % BASE);
            carry = prod / BASE;
        }
        res.trimLeadingZeros();
        return res;
    }

    template<typename D, typename DD>
    std::pair<BigInt<D, DD>, BigInt<D, DD>> BigInt<D, DD>::divideMagnitude(const BigInt& a, const BigInt& b)
    {
        if (b.isZero())
        {
            throw std::runtime_error("Division by zero");
        }
        if (compareMagnitude(a, b) < 0)
        {
            return {BigInt(0), a};
        }
        if (b == BigInt(1))
        {
            return {a, BigInt(0)};
        }
        if (b.digits_.size() == 1)
        {
            BigInt q = a;
            D r = divideShort(q, b.digits_[0]);
            return {q, BigInt(r)};
        }

        BigInt q, r;
        q.digits_.resize(a.digits_.size(), 0);

        for (int i = a.digits_.size() - 1; i >= 0; --i)
        {
            r.digits_.insert(r.digits_.begin(), 1, a.digits_[i]);
            r.trimLeadingZeros();

            if (compareMagnitude(r, b) < 0)
            {
                q.digits_[i] = 0;
                continue;
            }

            digit_type low = 1, high = BASE - 1, q_digit = 0;
            while (low <= high)
            {
                digit_type mid = low + (high - low) / 2;
                if (mid == 0)
                {
                    break;
                }
                BigInt temp_product = multiplyShort(b, mid);
                if (compareMagnitude(temp_product, r) <= 0)
                {
                    q_digit = mid;
                    low = mid + 1;
                }
                else
                {
                    high = mid - 1;
                }
            }

            if (q_digit > 0)
            {
                r -= multiplyShort(b, q_digit);
            }
            q.digits_[i] = q_digit;
        }

        q.trimLeadingZeros();
        r.trimLeadingZeros();
        return {q, r};
    }

    template <typename D, typename DD>
    bool operator<(const BigInt<D, DD>& a, const BigInt<D, DD>& b)
    {
        a.assertEagerState();
        b.assertEagerState();

        if (a.is_negative_ != b.is_negative_)
        {
            return a.is_negative_;
        }

        int mag_cmp = BigInt<D, DD>::compareMagnitude(a, b);
        if (a.is_negative_)
        {
            return mag_cmp > 0;
        }
        return mag_cmp < 0;
    }

    template <typename D, typename DD>
    bool operator==(const BigInt<D, DD>& a, const BigInt<D, DD>& b)
    {
        a.assertEagerState();
        b.assertEagerState();

        if (a.isZero() && b.isZero())
        {
            return true;
        }
        return a.is_negative_ == b.is_negative_ && BigInt<D, DD>::compareMagnitude(a, b) == 0;
    }
}