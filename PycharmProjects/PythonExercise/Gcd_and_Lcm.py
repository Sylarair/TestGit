# _*_ coding: utf-8 _*_
#import math


def gcd(m, n):
	if m < n:
		x, y = m, n
	else:
		x, y = n, m

	while x != 0:
		x, y = y % x, x

	return y


def lcm(m, n):
	_gcd = gcd(m, n)
	_lcm = m * n / _gcd

	return _lcm


def main():
	return gcd(12, 16), lcm(12, 16)

#gcd(12, 16)
main()
# if __name__ == '__main__':
# 	main()
