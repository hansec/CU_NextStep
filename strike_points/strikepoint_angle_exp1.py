import matplotlib.pyplot as plt
import numpy as np
from gfile_reader import gfile_read
from calc_bdry import calc_bdry
from calc_bdry import calc_xpoints
from calc_bdry import calc_strikepts
from calc_bdry import calc_alphas


#Read in a dictionary of the gfile information
gfile_dict1 = gfile_read(file0='d_0.5', silent=True)

Psi_rad = gfile_dict1['PSIRZ']
Psi = (2.0*np.pi)*Psi_rad  ##Convert from Wb/Rad to Wb
R_grid = gfile_dict1['R_grid']
Z_grid = gfile_dict1['Z_grid']
Psi_bdry = gfile_dict1['SIBRY'] #In wb/Rad

R_maxis = gfile_dict1['RMAXIS']
Z_maxis = gfile_dict1['ZMAXIS']

Limitr_R = gfile_dict1['RLIM']
Limitr_Z = gfile_dict1['ZLIM']
nlim = int(gfile_dict1['LIMITR'])
Limitr_R = Limitr_R[0:nlim]
Limitr_Z = Limitr_Z[0:nlim]

#Calculate the dictionary with the plasma boundary and sepratrix coordinates 
bdry_dict1 = calc_bdry(np.linspace(-np.pi, np.pi, 361), R_maxis, Z_maxis, R_grid, Z_grid,\
Psi_rad, Psi_bdry, xpt2=False, show_calc=False, reducebdry=0.0)# reducebdry=1.0)

#Calculate the dictionary with the X-point coordinates
xpoint_dict1 = calc_xpoints(R_grid, Z_grid, Psi_rad, Psi_bdry, Psi_restrict=False, xpt2_dict=False, show_calc=False)


#Calculate the upper strike pointc coordinates
try:
    upstrkpt_dict1 = calc_strikepts(bdry_dict1['upper_sep'][0], bdry_dict1['upper_sep'][1], Limitr_R, Limitr_Z,\
	xpoint_dict1['upper_xpt'][0], xpoint_dict1['upper_xpt'][1], calc_angles=False, show_calc=False)

    x_strike_up1 = [float(upstrkpt_dict1['left_strkpt'][0]), float(upstrkpt_dict1['right_strkpt'][0])]
    y_strike_up1 = [float(upstrkpt_dict1['left_strkpt'][1]), float(upstrkpt_dict1['right_strkpt'][1])]

except:
	x_strike_up1 = [0.0, 0.0]
	y_strike_up1 = [0.0, 0.0]
	print('Error in Finding Upper Strike Points')

#Calculate the upper strike point angles (3D)
upstrkpt_alpha_dict1 = calc_alphas(bdry_dict1['upper_sep'][0], bdry_dict1['upper_sep'][1], Limitr_R, Limitr_Z,\
gfile_dict1, upstrkpt_dict1, n_phi=0.0, show_calc=False)

print(upstrkpt_alpha_dict1)

#Calculate the lower strike pointc coordinates
try:
	lowstrkpt_dict1 = calc_strikepts(bdry_dict1['lower_sep'][0], bdry_dict1['lower_sep'][1], Limitr_R, Limitr_Z,\
	xpoint_dict1['lower_xpt'][0], xpoint_dict1['lower_xpt'][1], calc_angles=True, show_calc=False)

	x_strike_low1 = [float(lowstrkpt_dict1['left_strkpt'][0]), float(lowstrkpt_dict1['right_strkpt'][0])]
	y_strike_low1 = [float(lowstrkpt_dict1['left_strkpt'][1]), float(lowstrkpt_dict1['right_strkpt'][1])]
except:
	x_strike_low1 = [0.0, 0.0]
	y_strike_low1 = [0.0, 0.0]
	print('Error in Finding Lower Strike Points')

#Calculate the lower strike point angles (3D)
lowstrkpt_alpha_dict1 = calc_alphas(bdry_dict1['lower_sep'][0], bdry_dict1['lower_sep'][1], Limitr_R, Limitr_Z,\
gfile_dict1, lowstrkpt_dict1, n_phi=0.0, show_calc=False)

print(lowstrkpt_alpha_dict1)

#Plotting 
#====================================================
plt.figure(figsize=(5,8))
plt.plot(Limitr_R, Limitr_Z, color='k') #Plot the limiter contour

plt.plot(bdry_dict1['plasma_bdry'][0], bdry_dict1['plasma_bdry'][1], color='tab:blue')
plt.plot(bdry_dict1['upper_sep'][0], bdry_dict1['upper_sep'][1], color='tab:orange')
plt.plot(bdry_dict1['lower_sep'][0], bdry_dict1['lower_sep'][1], color='tab:green')

#Plot the xpoints
plt.plot(xpoint_dict1['lower_xpt'][0], xpoint_dict1['lower_xpt'][1], marker='o', color='k', markersize=3)
plt.plot(xpoint_dict1['upper_xpt'][0], xpoint_dict1['upper_xpt'][1], marker='o', color='k', markersize=3)

#Plot the strike points
plt.plot(upstrkpt_dict1['left_strkpt'][0], upstrkpt_dict1['left_strkpt'][1], marker='o', color='tab:red', markersize=3)
plt.plot(upstrkpt_dict1['right_strkpt'][0], upstrkpt_dict1['right_strkpt'][1], marker='o', color='tab:red', markersize=3)

plt.plot(lowstrkpt_dict1['left_strkpt'][0], lowstrkpt_dict1['left_strkpt'][1], marker='o', color='tab:red', markersize=3)
plt.plot(lowstrkpt_dict1['right_strkpt'][0], lowstrkpt_dict1['right_strkpt'][1], marker='o', color='tab:red', markersize=3)

#Label the strike point incidence angles
plt.text(upstrkpt_dict1['left_strkpt'][0], upstrkpt_dict1['left_strkpt'][1]+0.001,\
str(round(upstrkpt_alpha_dict1['leftstrkpt_alpha'], 3))+' deg', color='r', fontsize=6)

plt.text(upstrkpt_dict1['right_strkpt'][0]+0.01, upstrkpt_dict1['right_strkpt'][1]-0.01,\
str(round(upstrkpt_alpha_dict1['rightstrkpt_alpha'], 3))+' deg', color='r', fontsize=6)

plt.text(lowstrkpt_dict1['left_strkpt'][0], lowstrkpt_dict1['left_strkpt'][1]-0.02,\
str(round(upstrkpt_alpha_dict1['leftstrkpt_alpha'], 3))+' deg', color='r', fontsize=6)

plt.text(lowstrkpt_dict1['right_strkpt'][0]+0.01, lowstrkpt_dict1['right_strkpt'][1]+0.01,\
str(round(lowstrkpt_alpha_dict1['rightstrkpt_alpha'], 3))+' deg', color='r', fontsize=6)


plt.plot(R_maxis, Z_maxis, marker='+', color='tab:blue')

plt.gca().set_aspect('equal')
plt.xlabel('R (m)')
plt.ylabel('Z (m)')
plt.show()