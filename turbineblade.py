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
from scipy import integrate

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

	# draws 2 characteristic lines. the angle starts from v1(compression wave) and v2(expansion wave).
	def draw_lines(self,v1,v2):
		
		v1 = math.radians(v1)
		v2 = math.radians(v2)

		counter = 1000
		x0,y0,x1,y1,R_x,R_y,R_x_min,R_y_min = [],[],[],[],[],[],[],[]
		for num in range(1,counter):
			Rstar = self.Rstar_min + num*1/counter 
			if Rstar > 1: 
				Rstar = 1
			theta1 = (self.chara_line(Rstar) - self.const) + v1
			theta2 = -(self.chara_line(Rstar) - self.const) + v2 
			x0 += [(self.chara_x(Rstar,theta1))]
			y0 += [(self.chara_y(Rstar,theta1))]
			x1 += [(self.chara_x(Rstar,theta2))]
			y1 += [(self.chara_y(Rstar,theta2))]

			if Rstar == 1: break

#		for num in range(0,360):
#			R_x += [-math.sin(math.radians(num))]
#			R_y += [math.cos(math.radians(num))]
#			R_x_min += [-math.sin(math.radians(num))*self.Rstar_min]
#			R_y_min += [math.cos(math.radians(num))*self.Rstar_min]

		#test map func. Exactly same function as above "for" loop
		num = 360
		R_x = list(map(lambda theta: -math.sin(math.radians(theta)),range(0,num)))
		R_y = list(map(lambda theta: math.cos(math.radians(theta)),range(0,num)))
		R_x_min = list(map(lambda theta: -math.sin(math.radians(theta))*self.Rstar_min,range(0,num)))
		R_y_min = list(map(lambda theta: math.cos(math.radians(theta))*self.Rstar_min,range(0,num)))

		plt.plot(x0,y0)
		plt.plot(x1,y1)
		plt.plot(R_x,R_y)
		plt.plot(R_x_min,R_y_min)
		plt.xlim(-1,1)
		plt.ylim(-1,1)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()

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
		plt.xlim(-1,1)
		plt.ylim(-1,1)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()

	# top arc angle is 0 
	# v1 must be smaller than v2
	def get_R(self,v1,v2):

		v1 = math.radians(v1)
		v2 = math.radians(v2)
		Rstar = 1	#initial Rstar
		param = []
		temp_x = 1
		temp_y = 1
		myu_check = math.radians(90)

		#this for is very stupid. intersecting angle is always half the delta
		for num in range(1,10): #decimal precision of Rstar
			while(1):
				theta1 = (self.chara_line(Rstar) - self.const) + v1
				theta2 = -(self.chara_line(Rstar) - self.const) + v2 
				x0 = self.chara_x(Rstar,theta1)
				y0 = self.chara_y(Rstar,theta1)
				x1 = self.chara_x(Rstar,theta2)
				y1 = self.chara_y(Rstar,theta2)
#				print(Rstar,math.degrees(theta1),math.degrees(theta2))
#				print(x0,y0,x1,y1,Rstar)
				if x1 - x0 < 0: 
					if(temp_x - x1 == 0):		#this if is to avoid my fucking bug
						myu_check = math.radians(90)		#no meaning for the value. 
					else:
						myu_check = math.atan((temp_y - y1)/(temp_x - x1))
#						if myu_check == 0:
#							print(theta1,theta2)
					temp_x = x1
					temp_y = y1
					Rstar = Rstar + 1/(10**num) 

					break

				Rstar = Rstar - 1/(10**num)

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

	def get_Ru(self,vu):
		return self.get_R(-vu,vu)[0]

	#In the papaer, theta are defines around 90-120deg
	def get_upper_arc(self,Rstar,theta_in,theta_out,vu,vin,vout):
		
		x,y = [],[]

		beta_in = -(theta_in-vu+vin)
		beta_out = theta_out-vu+vout

		theta = np.arange(beta_in,beta_out,0.1)
		x = list(map(lambda theta: Rstar*math.sin(math.radians(theta)),theta))
		y = list(map(lambda theta: Rstar*math.cos(math.radians(theta)),theta))
		
		print(beta_in,beta_out)
		print(x[0],y[0])
		print(x[-1],y[-1])
		plt.plot(x,y)
		plt.xlim(-2,2)
		plt.ylim(-2,2)
		plt.gca().set_aspect('equal', adjustable='box')
#		plt.show()

	def get_lower_arc(self,Rstar,theta_in,theta_out,vl,vin,vout,shift = 0):
		
		x,y = [],[]

		alpha_in = -(theta_in-vin+vl)
		alpha_out = theta_out-vout+vl

		theta = np.arange(alpha_in,alpha_out,0.1)
		x = list(map(lambda theta: Rstar*math.sin(math.radians(theta)),theta))
		y = list(map(lambda theta: shift + Rstar*math.cos(math.radians(theta)),theta))
		
		print(alpha_in,alpha_out)
		print(x[0],y[0])
		print(x[-1],y[-1])
		
		plt.plot(x,y)
		plt.xlim(-2,2)
		plt.ylim(-2,2)
		plt.gca().set_aspect('equal', adjustable='box')
#		plt.show()


	# defines upper convex arc coordinates
	# ve is the entry angle
	# vu is the value to define the upper convex radius
	def upper_convex(self,vu,ve,theta):
		
		x,y = [],[]
		Xstar,Ystar = [],[]
		R_x,R_y = [],[]
		R_x_min,R_y_min = [],[]

		xtmp = 0
		ytmp = self.get_Ru(vu)

		for num in range(0,ve*2):
			Xstar_b = xtmp
			Ystar_b = ytmp
			Rstar,Xstar_a,Ystar_a,myu_check = self.get_R(-vu,vu-num)
			myu = self.get_myu(self.get_mach(self.get_Mstar(Rstar)))
			a1 = math.tan(myu+math.radians(num/2))
#			a1 = math.tan(myu_check)
			b1 = Ystar_a - a1 * Xstar_a
			a2 = math.tan(math.radians(num/2))
#			a2 = math.tan(myu_check-myu)
			b2 = Ystar_b - a2 * Xstar_b
			
#			print(Rstar,Xstar_a,Ystar_a,math.degrees(myu_check),math.degrees(myu))
			
			xtmp = ((b2 - b1) / (a1 - a2))
			ytmp = xtmp * a2 + b2

#			print(num/2,Rstar,math.degrees(myu_check),math.degrees(myu),xtmp,ytmp)

			rotx,roty = self.rotate(xtmp,ytmp,theta)

			x += [(rotx)]
			y += [(roty)]
			Xstar += [(Xstar_a)]
			Ystar += [(Ystar_a)]

			xcomp = np.arange(-1,0,0.001)
			ycomp = a1*xcomp + b1

#			plt.plot(xcomp,ycomp)

		for num in range(0,90):
			R_x += [-math.sin(math.radians(num))]
			R_y += [math.cos(math.radians(num))]
			R_x_min += [-math.sin(math.radians(num))*self.get_Ru(vu)]
			R_y_min += [math.cos(math.radians(num))*self.get_Ru(vu)]


		plt.plot(x,y)
#		plt.plot(Xstar,Ystar)
#		plt.plot(R_x,R_y)
#		plt.plot(R_x_min,R_y_min)
#		plt.xlim(-1,0)
#		plt.ylim(0,1)
#		plt.gca().set_aspect('equal', adjustable='box')
#		plt.show()

		return x[-1],y[-1]

	# defines lower concave arc coordinates
	# Precaution. there are still calculation uncertainties.
	# ve is the entry angle
	def lower_concave(self,v1,ve,theta,shift = 0):
		
		xinit = -math.sin(math.radians(v1))		
		yinit= math.cos(math.radians(v1))		

		xtmp = xinit
		ytmp = yinit

		x,y = [],[]
		Xstar,Ystar = [],[]
		R_x,R_y = [],[]

		#v1 is the origin angle. ve is the end angle
		for num in range(v1,ve*2):	

			Xstar_b = xtmp
			Ystar_b = ytmp
			Rstar,Xstar_a,Ystar_a,myu_check = self.get_R(-num,v1)
			myu = self.get_myu(self.get_mach(self.get_Mstar(Rstar)))
			a1 = math.tan(-myu+math.radians(num/2))
			b1 = Ystar_a - a1 * Xstar_a
			a2 = math.tan(math.radians(num/2))
			b2 = Ystar_b - a2 * Xstar_b

			xtmp = ((b2 - b1) / (a1 - a2))
			ytmp = xtmp * a2 + b2

			print(num/2,Rstar,math.degrees(myu_check),math.degrees(myu),xtmp,ytmp)
			
			rotx,roty = self.rotate(xtmp,ytmp,theta)

			x += [(rotx)]
			y += [(roty + shift)]
			Xstar += [(Xstar_a)]
			Ystar += [(Ystar_a)]

			xcomp = np.arange(-1,0,0.001)
			ycomp = a1*xcomp + b1

#			plt.plot(xcomp,ycomp)

		for num in range(0,90):
			R_x += [-math.sin(math.radians(num))]
			R_y += [math.cos(math.radians(num))]

		plt.plot(x,y)
#		plt.plot(Xstar,Ystar)
#		plt.plot(R_x,R_y)
#		plt.xlim(-1,0)
#		plt.ylim(0,1)
#		plt.gca().set_aspect('equal', adjustable='box')
#		plt.show()

		return x[-1],y[-1]

	def get_Q(self,Rlstar,Rustar):
		
		Mlstar = 1/Rlstar
		Mustar = 1/Rustar

		Mach = lambda Mach:(((self.gamma+1) / 2 - (self.gamma-1) * Mach**2 / 2)**(1 / (self.gamma - 1)))/Mach

		Q1 = Mlstar * Mustar / (Mustar - Mlstar)
		Q2 = integrate.quad(Mach,Mlstar,Mustar)

		Q = Q1*Q2[0]

		return Q

	def get_Gstar(self,vl,vu,theta_in):

		Ae = 0
		Rlstar = self.get_Ru(vl)
		Rustar = self.get_Ru(vu)
		Q = self.get_Q(Rlstar,Rustar)
		
		Astar = (Rlstar - Rustar) / Q

		print(Astar)

#		Gstar = Ae/Astar * Q * (Rl - Ru) / math.cos(theta_in)
		print("daradara")


	def rotate(self,x,y,angle):

		theta = math.radians(angle)

		a = np.array((
			(math.cos(theta),-math.sin(theta)),
			(math.sin(theta),math.cos(theta))
			))

		b = np.array((x,y))

		return np.dot(a,b)

	def straight_line(self,theta,x,y,targetx):

		targety = math.tan(math.radians(theta)) * targetx + y - math.tan(math.radians(theta)) * x

		return targety


if __name__ == "__main__":

	gamma = 1.4
	vin = 30
	vout = 20
	vl = 0
	vu = 30
	theta_in_upper = 30
	theta_in_lower = 60
	theta_out_upper = 60
	theta_out_lower = 60

	print("Design Supersonic Turbine")

	f = Blade(gamma)

#	f.get_upper_arc(f.get_Ru(vu),theta_in_upper,theta_out_upper,vu,vin,vout)
#	f.get_lower_arc(f.get_Ru(vl),theta_in_lower,theta_out_lower,vl,vin,vout)
#	plt.show()
	
#	a = f.rotate(1,0,45)
#	f.get_Gstar(vl,vu,theta_in_upper)

	xlow,ylow = f.lower_concave(vl,vin,theta_in_lower-vin)
	xup,yup = f.upper_convex(vu,vin,theta_in_upper)
	f.get_lower_arc(f.get_Ru(vl),theta_in_lower,theta_out_lower,vl,vin,vout)
	f.get_upper_arc(f.get_Ru(vu),theta_in_upper,theta_out_upper,vu,vin,vout)
	newy = f.straight_line(90-vin,xup,yup,xlow)
	plt.plot([xlow,xup],[newy,yup])
	shift = -abs(ylow - newy)
	f.lower_concave(vl,vin,theta_in_lower-vin,shift)
	f.get_lower_arc(f.get_Ru(vl),theta_in_lower,theta_out_lower,vl,vin,vout,shift)
	plt.show()
#	f.draw_lines(-10,20)
#	f.get_R(-10,10)
#	f.lower_concave(0,30)
#	f.upper_convex(40,35)
#	f.plot_chara(360,5)
	print("finish")
