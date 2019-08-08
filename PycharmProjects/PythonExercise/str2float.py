# _*_ coding: utf-8 _*_
import sys


def add(x, y):
	#if isinstance(_input, str):
	return 10 * x + y

def str2float(_input):
	str2int_dict = {}
	for i in range(0, 10):
		str2int_dict[str(i)] = int(i)

	return str2int_dict[_input]
	#str_list = _input.split('')
	#point_pos = len(str_list) - str_list.index('.') - 1

	#return

def get_result(_string):
	if isinstance(_string, str):
		point_pos = _string.index('.')
		divid_num = 10 ** (len(_string) - point_pos -1)
		print point_pos
		__string = _string[:point_pos] + _string[point_pos + 1:]
		num_1 = reduce(add, map(str2float, __string))
		print num_1
		num_2 = float(num_1) / divid_num
		print num_2

		return num_2
	else:
		print 'Error! Please input correct string!'
		sys.exit(1)

def main():
	get_result('1234.5')