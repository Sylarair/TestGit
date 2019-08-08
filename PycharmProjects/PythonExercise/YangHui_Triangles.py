# _*_ coding : utf-8 _*_
import sys
def triangles(n):
	L = [1]
	# L_new = [1,1]
	if n < 1:
		print 'Error! Please check your input number!'
		sys.exit(1)
	elif n == 1:
		print L
	elif n == 2:
		print L
		print [1, 1]
	else:
		for i in range(0, n):
			if i == 0:
				print L
			elif i == 1:
				L = [1, 1]
				print L
			else:
				# for j in range(0,i-1):
				L = [1] + [L[j - 1] + L[j] for j in range(1, len(L))] + [1]
				print L

def main():
	triangles(20)

if __name__ == '__main__':
	main()