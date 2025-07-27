//
// Created by Ext Culith on 2025/7/27.
//

#include "../include/BigInt.h"
#include <string>

BigInt::BigInt(std::string str) : isMinus(str[0] == '-')
{
	vec.reserve(str.size() + 20);
	int i = 0;
	if (isMinus) i = 1;
	for (int j = i; j < str.size(); ++j)
	{
		if (!isdigit(str[j]))
			throw bad_input_string();
	}
	for (auto it1 = str.rbegin(), it2 = str.rend() - i; it1 != it2; ++it1)  //may be a joke
		vec.push_back(*it1 - '0');
	syncSize();
}
void BigInt::setStr(std::string str)
{
	vec.resize(0);
	vec.reserve(str.size() + 20);
	int i = 0;
	if (str[0] == '-') i = 1, isMinus = true;
	for (int j = i; j < str.size(); ++j)
	{
		if (!isdigit(str[j]))
			throw bad_input_string();
	}
	for (auto it1 = str.rbegin(), it2 = str.rend() - i; it1 != it2; ++it1)
		vec.push_back(*it1 - '0');
	syncSize();
}
void BigInt::syncSize()
{
	for (int i = vec.size() - 1;; --i)
	{
		if (vec[i] != 0) break;
		vec.pop_back();
		if (vec.size() == 1 || vec.empty())
			return;
	}
}
void BigInt::inverse()
{
	for (auto& x : vec)
		x *= -1;
}
void BigInt::syncMinus()  //More efficient within flattened vector
{
	syncSize();
	vec.push_back(0);
	for (int i = 0; i < vec.size() - 1; ++i)
	{
		if (vec[i] < 0)
		{
			vec[i] += 10;
			vec[i + 1] -= 1;
		}
	}
	if (vec[vec.size() - 1] >= 0);
	else
	{
		isMinus = !isMinus;
		for (int i = 0; i < vec.size() - 1; ++i)
		{
			if (vec[i] != 0)
			{
				vec[i] = 10 - vec[i];
				vec[i + 1] += 1;
			}
		}
	}
	syncSize();
}
BigInt::vecType BigInt::getInverse() const
{
	vecType tmp;
	for (auto& x : vec)
		tmp.push_back(x * -1);
	return tmp;
}
void BigInt::sync()
{
	int tmp;
	syncSize();
	vec.push_back(0);
	for (int i = 0; i < vec.size() - 1; ++i)
	{
		if (vec[i] < 0)
		{
			tmp = -vec[i];
			if (tmp % 10 == 0)
			{
				vec[i] = 0;
				vec[i + 1] -= tmp / 10;
			}
			else
			{
				vec[i] = 10 - tmp % 10;
				vec[i + 1] -= tmp / 10 + 1;
			}
		}
		else
		{
			vec[i + 1] += vec[i] / 10;
			vec[i] %= 10;
		}
	}
	tmp = vec[vec.size() - 1];
	if (tmp >= 0)
	{
		int n = 0;
		for (short i = vec[vec.size() - 1]; i > 0; ++n) i /= 10;
		vec.resize(vec.size() + n - 1, 0);
		for (int k = 0; tmp != 0; ++k)
		{
			vec[vec.size() - n + k] = tmp % 10;
			tmp /= 10;
		}
	}
	else
	{
		isMinus = !isMinus;
		tmp *= -1;
		int n = 0;
		for (short i = tmp; i > 0; ++n) i /= 10;
		vec.resize(vec.size() + n - 1, 0);
		for (int k = 0; tmp != 0; ++k)
		{
			vec[vec.size() - n + k] = -(tmp % 10);
			tmp /= 10;
		}
		if (vec[vec.size() - 1] != -1)
		{
			vec.push_back(0);
			for (int k = n + 1; k > 1; --k)
			{
				vec[vec.size() - k] += 10;
				vec[vec.size() - k + 1] -= 1;
			}
		}
		for (int i = 0; i < vec.size() - 1; ++i)
		{
			if (vec[i] != 0)
			{
				vec[i] = 10 - vec[i];
				vec[i + 1] += 1;
			}
		}
	}
	syncSize();
}
BigInt& BigInt::addNoFlat(const BigInt& num)
{
	if (isMinus) inverse();
	if (vec.size() < num.vec.size())
		vec.resize(num.vec.size());
	int i = 0;
	for (auto& x : num.vec)
	{
		vec[i] += x;
		++i;
	}
	return *this;
}
BigInt& BigInt::minusNoFlat(const BigInt& num)
{
	if (isMinus) inverse();
	if (vec.size() < num.vec.size())
		vec.resize(num.vec.size());
	int i = 0;
	for (auto& x : num.vec)
	{
		vec[i] -= x;
		++i;
	}
	return *this;
}
BigInt& BigInt::addNoSign(const BigInt& num)
{
	if (vec.size() <= num.vec.size())
		vec.resize(num.vec.size() + 2, 0);
	int i = 0;
	short tmp = 0;
	for (auto& x : num.vec)
	{
		tmp = vec[i] + x;
		if (tmp >= 10)
		{
			vec[i] = tmp % 10;
			vec[i + 1] += 1;
			for (int k = 1; vec[i + k] == 10; ++k)
			{
				vec[i + k + 1] += 1;
				vec[i + k] = 0;
			}
		}
		else vec[i] = tmp;
		++i;
	}
	syncSize();
	return *this;
}
BigInt& BigInt::minusNoSign(const BigInt& num)
{
	auto tmp = num.getInverse();
	if (vec.size() <= num.vec.size())
		vec.resize(num.vec.size() + 1, 0);
	int i = 0;
	for (auto& x : tmp)
	{
		vec[i] += x;
		++i;
	}
	syncMinus();
	return *this;
}
BigInt& BigInt::multipleNoSign(const BigInt& num)
{
	int n = 0;
	vecType tmp = vec;
	vec.resize(num.vec.size() + vec.size(), 0);
	int k = 0;
	for (auto& x : num.vec)
	{
		switch (n)
		{
		case 0:
		{
			k = 0;
			for (auto& x1 : tmp)
			{
				vec[k] = x1 * x;
				++k;
			}
			break;
		}
		default:
		{
			k = n;
			for (auto& x1 : tmp)
			{
				vec[k] += x1 * x;
				++k;
			}
		}
		}
		++n;
	}
	sync();
	return *this;
}
BigInt& BigInt::multipleInt(int num)
{
	auto mI = [this](const int& num) -> void
		{
			for (auto& x : vec)
				x *= num;
			sync();
		};
	if (isMinus && num < 0)
	{
		num *= -1;
		mI(num);
		isMinus = !isMinus;
	}
	if (!isMinus && num < 0)
	{
		num *= -1;
		mI(num);
		isMinus = !isMinus;
	}
	if (isMinus && num >= 0)
		mI(num);
	else
		mI(num);
	return *this;
}
BigInt& BigInt::operator+=(const BigInt& num)
{
	if (isMinus && num.isMinus)
		addNoSign(num);
	if (!isMinus && num.isMinus)
		minusNoSign(num);
	if (isMinus && !num.isMinus)
	{
		minusNoSign(num);
		isMinus = !isMinus;
	}
	else
		addNoSign(num);
	return *this;
}
BigInt& BigInt::operator-=(const BigInt& num)
{
	if (isMinus && num.isMinus)
	{
		minusNoSign(num);
		isMinus = !isMinus;
	}
	if (!isMinus && num.isMinus)
		addNoSign(num);
	if (isMinus && !num.isMinus)
		addNoSign(num);
	else
		minusNoSign(num);
	return *this;
}
bool BigInt::lessNoSign(const BigInt& num) const
{
	if (vec.size() > num.vec.size())
		return false;
	if (vec.size() < num.vec.size())
		return true;
	for (int n = 0; n < vec.size(); ++n)
	{
		if (vec[n] < num.vec[n])
			return true;
		if (vec[n] > num.vec[n])
			return false;
	}
	return false;
}
bool BigInt::notMoreNoSign(const BigInt& num) const
{
	if (vec.size() > num.vec.size())
		return false;
	if (vec.size() < num.vec.size())
		return true;
	for (int n = 0; n < vec.size(); ++n)
	{
		if (vec[n] < num.vec[n])
			return true;
		if (vec[n] > num.vec[n])
			return false;
	}
	return true;
}
bool BigInt::operator<(const BigInt& num)
{
	if (isMinus && num.isMinus)
		return !lessNoSign(num);
	if (!isMinus && num.isMinus)
		return false;
	if (isMinus && !num.isMinus)
		return true;
	else
		return lessNoSign(num);
}
BigInt& BigInt::operator*=(const BigInt& num)
{
	if (isMinus && num.isMinus)
	{
		multipleNoSign(num);
		isMinus = !isMinus;
	}
	if (!isMinus && num.isMinus)
	{
		multipleNoSign(num);
		isMinus = !isMinus;
	}
	if (isMinus && !num.isMinus)
		multipleNoSign(num);
	else
		multipleNoSign(num);
	return *this;
}
BigInt& BigInt::divideNoSign(const BigInt& num)
{
	if (num.vec[num.vec.size() - 1] == 1)
	{
		bool status = true;
		for (int n = vec.size() - 2; n >= 0; --n)
			if (vec[n] != 0) status = false;
		if (status == true)
		{
			vec.erase(vec.begin());
			return *this;
		}
	}
	long len = num.vec.size();
	BigInt tmp = 0;
	tmp.vec.resize(vec.size() - num.vec.size() - 1, 0);
	for (auto& x : num.vec)
	{
		tmp.vec.push_back(x);
	}
	vecType result(vec.size() - len, 0);
	for (long r = vec.size() - num.vec.size(); r > 0; --r)
	{
		while (tmp.notMoreNoSign(*this))
		{
			this->minusNoSign(tmp);
			result[r - 1] += 1;
		}
		tmp.divideNoSign(10);
		syncSize();
	}
	vec = result;
	return *this;
}
BigInt& BigInt::operator/=(const BigInt& num)
{
	if (isMinus && num.isMinus)
	{
		divideNoSign(num);
		isMinus = !isMinus;
	}
	if (!isMinus && num.isMinus)
	{
		divideNoSign(num);
		isMinus = !isMinus;
	}
	if (isMinus && !num.isMinus)
		divideNoSign(num);
	else
		divideNoSign(num);
	return *this;
}
BigInt& BigInt::modNoSign(const BigInt& num)
{
	if (num.vec[num.vec.size() - 1] == 1)
	{
		bool status = true;
		for (int n = vec.size() - 2; n >= 0; --n)
			if (vec[n] != 0) status = false;
		if (status == true)
		{
			vec.erase(vec.begin());
			return *this;
		}
	}
	long len = num.vec.size();
	BigInt tmp = 0;
	tmp.vec.resize(vec.size() - num.vec.size() - 1, 0);
	for (auto& x : num.vec)
	{
		tmp.vec.push_back(x);
	}
	for (long r = vec.size() - num.vec.size(); r > 0; --r)
	{
		while (tmp.notMoreNoSign(*this))
		{
			this->minusNoSign(tmp);
		}
		tmp.divideNoSign(10);
		syncSize();
	}
	return *this;
}
std::string BigInt::ostr() const
{
	std::string str;
	if (isMinus) str.push_back('-');
	for (auto it1 = vec.rbegin(), it2 = vec.rend(); it1 != it2; ++it1)
		str.push_back(*it1 + '0');
	return str;
}