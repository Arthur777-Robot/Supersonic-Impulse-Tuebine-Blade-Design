# -*- coding: utf-8 -*-
'''
The MIT License (MIT)
Copyright (c) 2017 Seiji Arther Murakami

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

plt.close('all')


# Characteristic line function

class Blade():
	def __init__(self,gamma):
		print('Initializing parameter...')
		self.gamma = gamma
		self.Rstar_min = math.sqrt((self.gamma - 1)/(self.gamma + 1))
		self.const = self.chara_line(1)

	def chara_line(self,Rstar):

		fai = 0.5 * (math.sqrt((self.gamma + 1)/(self.gamma - 1)) * math.asin((self.gamma - 1) / Rstar**2 - self.gamma) + math.asin((self.gamma + 1) *  Rstar**2 - self.gamma))

		return fai

	def chara_x(self,Rstar,theta):
		return Rstar*math.sin(theta) 

	def chara_y(self,Rstar,theta):
		return Rstar*math.cos(theta) 

	def plot_chara(self,angle,step):

		x1,y1,x2,y2 = [],[],[],[]

		counter = 1000
		for j in range(0,angle,step):
			del x1[:],y1[:],x2[:],y2[:]
			for i in range(counter):
				Rstar = self.Rstar_min + i *1/counter 
				if Rstar > 1: 
					Rstar = 1
				fai = self.chara_line(Rstar)
				x1.append(self.chara_x(Rstar,fai - self.const + math.radians(j)))
				y1.append(self.chara_y(Rstar,fai - self.const + math.radians(j)))
				x2.append(self.chara_x(Rstar,-(fai - self.const - math.radians(j))))
				y2.append(self.chara_y(Rstar,-(fai - self.const - math.radians(j))))

				if Rstar == 1: break

			plt.plot(x1,y1,"r")
			plt.plot(x2,y2,"k")
		plt.show()

	def get_R(self,org,nyu):
		
		delta = math.radians(nyu)*2
		org = math.radians(org)
		Rstar = self.Rstar_min
		param = []
		temp_x = 1
		temp_y = 1

		#this for is very stupid. intersecting angle is always half the delta
		for num in range(1,10): #decimal precision of Rstar
			while(1):
				if Rstar >= 1:
					Rstar = 1	
				theta1 = self.chara_line(Rstar) - self.const + org
				theta2 = theta1 - delta
				x0 = self.chara_x(Rstar,-theta1)
				y0 = self.chara_y(Rstar,-theta1)
				x1 = self.chara_x(Rstar,theta2)
				y1 = self.chara_y(Rstar,theta2)
				if x0 - x1 > 0: 
					if(temp_x - x1 == 0):		#this if is to avoid my fucking bug
						myu_check = math.radians(90)		#no meaning for the value. 
					else:
						myu_check = math.atan((temp_y - y1)/(temp_x - x1))
					temp_x = x1
					temp_y = y1
					Rstar = Rstar - 1/(10**num) 

					break

				Rstar = Rstar + 1/(10**num)

		param.append(Rstar)
		param.append(x0)
		param.append(y0)
		param.append(myu_check)
		return param

	def get_myu(self,M):
		return math.asin(1/M)
	
	def get_Mstar(self,Rstar):
		return 1/Rstar

	def get_mach(self,Mstar):
		return math.sqrt((2 * Mstar**2)/((self.gamma + 1) - (self.gamma - 1) * Mstar**2))

	# v1 is the start angle,ve is end angle
	def lower_concave(self,v1,ve):

		myu = math.radians(2)
		fai = math.radians(4)
		Rstar,Xstar_a,Ystar_a = self.get_R(0,2)
		
		b1 = math.tan(myu + fai)
		c1 = -(Ystar_a + Xstar_a * math.tan(myu + fai))
		b2 = -math.tan(-fai)

		delta3 = b2 - b1
		
		f = delta3 + math.tan(-fai)
		g = delta3 + b1

		print(b1,c1,b2,delta3,f,g)

#		print(((c1 / f) - (c1 * b2) / (f * g)) ,  (b1 * math.tan(-fai)) , (f * g))
#		Xstar_b = ((c1 / f) - (c1 * b2) / (f * g)) / (1 - (b1 * math.tan(-fai)) / (f * g))
		Xstar_b = 0
		Ystar_b = ((b1 * Xstar_b * math.tan(-fai)) - c1 * b1) / g

		print(Xstar_b,Ystar_b)

		
	# vu is the start angle, Ve is the end angle
	def upper_convex(self,vu,ve):

		myu = math.radians(2)
		fai = math.radians(4)
		
		Rstar,Xstar_a,Ystar_a = self.get_R(60,2)

		b1 = -math.tan(myu-fai)
		c1 = Xstar_a * math.tan(myu - fai) - Ystar_a
		b2 = -math.tan(-fai)

		delta3 = b2 - b1

		f = delta3 + math.tan(-fai)
		g = delta3 - b2

		print(b1,c1,b2,delta3,f,g)

		Xstar_b = ((b1 * c1) / (f * g) + (c1 / f)) / (1.0 + (b2 * math.tan(-fai)) / (f * g))
		Ystar_b = ((b1 * c1) - b2 * Xstar_b * math.tan(-fai)) / g

		print(Xstar_b,Ystar_b)

	# defines lower concave arc coordinates
	# Precaution. there are still calculation uncertainties.
	def new_concave(self,fai):
		
		v1 = 0
		ve = 38

		fai = math.radians(fai)

		xtmp = 0		#initial Xstar_b = 0
		ytmp = 1		#initial Ystar_b = 1

		x,y = [],[]
		Xstar,Ystar = [],[]
		R_x,R_y = [],[]

		#v1 is the origin angle. ve is the end angle
		for num in range(v1 + 1,ve):	

			Xstar_b = xtmp
			Ystar_b = ytmp
			
			Rstar,Xstar_a,Ystar_a,myu_check = self.get_R(v1,num/2)
			myu = self.get_myu(self.get_mach(self.get_Mstar(Rstar)))
			a1 = math.tan(myu_check)
			b1 = Ystar_a - a1 * Xstar_a
			a2 = math.tan(math.radians(num))
			b2 = Ystar_b - a2 * Xstar_b

#			print(Ystar_a,a1*Xstar_a)
#			print(a1,b1,a2,b2,math.degrees(myu_check))
#			print(Xstar_a,Ystar_a)

#			print(Ystar_a)
#			print(math.degrees(math.atan(a1)),b1,math.degrees(math.atan(a2)),b2)

			xtmp = ((b2 - b1) / (a1 - a2))
			ytmp = xtmp * a2 + b2

#			print(Xstar_a,Ystar_a,math.degrees(myu_check))
			print(xtmp,ytmp)
			x += [float(xtmp)]
			y += [float(ytmp)]
			Xstar += [(Xstar_a)]
			Ystar += [(Ystar_a)]
			R_x += [-math.sin(math.radians(num))]
			R_y += [math.cos(math.radians(num))]

		plt.plot(x,y)
		plt.plot(Xstar,Ystar)
		plt.plot(R_x,R_y)
		plt.xlim(-1,0)
		plt.ylim(0,1)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()

if __name__ == "__main__":
	print("Design Supersonic Turbine")
	f = Blade(1.4)
	print("lower_concave")

	f.new_concave(1)
#	f.plot_chara(360,5)
	print("finish")
