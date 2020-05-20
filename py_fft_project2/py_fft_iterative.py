import math

def complex_dft(xr, xi, n):
	pi = 3.141592653589793
	rex = [0] * n
	imx = [0] * n
	for k in range(0, n):  # exclude n
		rex[k] = 0
		imx[k] = 0
	for k in range(0, n):  # for each value in freq domain
		for i in range(0, n):  # correlate with the complex sinusoid
			sr =  math.cos(2 * pi * k * i / n)
			si = -math.sin(2 * pi * k * i / n)
			rex[k] += xr[i] * sr - xi[i] * si
			imx[k] += xr[i] * si + xi[i] * sr
	return rex, imx

# FFT version based on the original BASIC program
def fft_basic(rex, imx, n):
	pi = 3.141592653589793
	m = int(math.log(n, 2))  # float to int
	j = n // 2

	# bit reversal sorting
	for i in range(1, n - 1):  # [1,n-2]
		if i >= j:
			# swap i with j
			rex[i], rex[j] = rex[j], rex[i]
			imx[i], imx[j] = imx[j], imx[i]
		k = n // 2
		while (1):
			if k > j:
				break
			j -= k
			k //= 2
		j += k

	for l in range(1, m + 1):  # each stage
		le = int(math.pow(2, l))  # 2^l
		le2 = le // 2
		ur = 1
		ui = 0
		sr =  math.cos(pi / le2)
		si = -math.sin(pi / le2)
		for j in range(1, le2 + 1):  # [1, le2] sub DFT
			for i in range(j - 1, n - 1, le):  #  for butterfly
				ip = i + le2
				tr = rex[ip] * ur - imx[ip] * ui
				ti = rex[ip] * ui + imx[ip] * ur
				rex[ip] = rex[i] - tr
				imx[ip] = imx[i] - ti
				rex[i] += tr
				imx[i] += ti
			tr = ur
			ur = tr * sr - ui * si
			ui = tr * si + ui * sr
