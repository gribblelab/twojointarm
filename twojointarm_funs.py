import numpy as np

arm_params = {
	'l': np.array([0.3384, 0.4554]),  # link lengths s,e (metres)
	'r': np.array([0.1692, 0.2277]),  # radius of gyration (metres)
	'm': np.array([2.10, 1.65]),      # mass (kg)
	'i': np.array([0.025, 0.075]),    # moment of inertia (kg*m*m)
}

def minjerk(H1,H2,t,n):
	# Given hand initial position H1=(x1,y1), final position H2=(x2,y2) and movement duration t,
	# and the number of points n,
	# Calculates the hand path H over time T that satisfies minimum-jerk.
	# Also returns derivatives Hd and Hdd.
	# From the paper:
	#   Flash T. & Hogan N. (1985) "The coordination of arm movements: an experimentally confirmed
	#        mathematical model" Journal of Neuroscience 5(7): 1688-1703.
	T = np.linspace(0,t,n)
	H = np.zeros((n,2))
	Hd = np.zeros((n,2))
	Hdd = np.zeros((n,2))
	t3,t4,t5 = t**3, t**4, t**5
	T2,T3,T4 = T**2, T**3, T**4
	Tt3,Tt4,Tt5 = (T/t)**3, (T/t)**4, (T/t)**5
	H[:,0]   =  H1[0] + (H1[0]-H2[0])*(15*Tt4 - 6*Tt5 - 10*Tt3)
	H[:,1]   =  H1[1] + (H1[1]-H2[1])*(15*Tt4 - 6*Tt5 - 10*Tt3)
	Hd[:,0]  = (H1[0] - H2[0])*(-30*T4/t5 + 60*T3/t4 - 30*T2/t3)
	Hd[:,1]  = (H1[1] - H2[1])*(-30*T4/t5 + 60*T3/t4 - 30*T2/t3)
	Hdd[:,0] = (H1[0] - H2[0])*(-120*T3/t5 + 180*T2/t4 - 60*T/t3)
	Hdd[:,1] = (H1[1] - H2[1])*(-120*T3/t5 + 180*T2/t4 - 60*T/t3)
	return (T,H,Hd,Hdd)

def joints_to_hand(A, arm_params):
	l1 = arm_params['l'][0]
	l2 = arm_params['l'][1]
	if A.ndim == 1:
		E = np.array([l1 * np.cos(A[0]), l1 * np.sin(A[0])])
		H = E + np.array([l2 * np.cos(A[0]+A[1]), l2 * np.sin(A[0]+A[1])])
	else:
		E = np.stack(([l1 * np.cos(A[:,0]), l1 * np.sin(A[:,0])]), axis=-1)
		H = E + np.stack(([l2 * np.cos(A[:,0]+A[:,1]), l2 * np.sin(A[:,0]+A[:,1])]), axis=-1)
	return (H,E)

def hand_to_joints(H, arm_params):
	l1 = arm_params['l'][0]
	l2 = arm_params['l'][1]
	ndim = np.ndim(H)
	if ndim == 1: # a single s,e joint angle pair
		A = np.zeros(2)
		A[1] = np.arccos((H[0]**2 + H[1]**2 - l1**2 - l2**2) / (2*l1*l2))
		A[0] = np.arctan2(H[1],H[0]) - np.arctan2(l2*np.sin(A[1]),l1+(l2*np.cos(A[1])))
	elif ndim == 2:  # time-varying s,e joint angles (nx2)
		n = np.shape(H)[0]
		A = np.zeros((n,2))
		A[:,1] = np.arccos((H[:,0]**2 + H[:,1]**2 - l1**2 - l2**2) / (2*l1*l2))
		A[:,0] = np.arctan2(H[:,1],H[:,0]) - np.arctan2(l2*np.sin(A[:,1]),l1+(l2*np.cos(A[:,1])))
	elif ndim == 3: # time varying (s,e) angles,vels,accels (nx2),(nx2),(nx2)
		H,Hd,Hdd = H
		n = np.shape(H)[0]
		A,Ad,Add = np.zeros((n,2)), np.zeros((n,2)), np.zeros((n,2))
		A[:,1] = np.arccos((H[:,0]**2 + H[:,1]**2 - l1**2 - l2**2) / (2*l1*l2))
		A[:,0] = np.arctan2(H[:,1],H[:,0]) - np.arctan2(l2*np.sin(A[:,1]),l1+(l2*np.cos(A[:,1])))
		for i in range(n):
		   J = arm_jacobian(A[i,:],arm_params)
		   Ad[i,:] = (np.linalg.inv(J).dot(Hd[i,:].T).T)
		   Jd = arm_jacobian_d(A[i,:],Ad[i,:],arm_params)
		   b = Hdd[i,:] - (Jd.dot(Ad[i,:].T)).T
		   Add[i,:] = (np.linalg.inv(J).dot(b.T)).T
		A = (A,Ad,Add)
	return A

def arm_jacobian(A, arm_params):
	l1 = arm_params['l'][0]
	l2 = arm_params['l'][1]
	J = np.zeros((2,2))
	J[0,0] = -l1*np.sin(A[0]) - l2*np.sin(A[0]+A[1])
	J[0,1] = -l2*np.sin(A[0]+A[1])
	J[1,0] = l1*np.cos(A[0]) + l2*np.cos(A[0]+A[1])
	J[1,1] = l2*np.cos(A[0]+A[1])
	return J

def arm_jacobian_d(A, Ad, arm_params):
	l1 = arm_params['l'][0]
	l2 = arm_params['l'][1]
	Jd = np.zeros((2,2))
	Jd[0,0] = -l1*np.cos(A[0])*Ad[0] - l2*(Ad[0] + Ad[1])*np.cos(A[0] + A[1])
	Jd[0,1] = -l2*(Ad[0] + Ad[1])*np.cos(A[0] + A[1])
	Jd[1,0] = -l1*np.sin(A[0])*Ad[0]  - l2*(Ad[0] + Ad[1])*np.sin(A[0] + A[1])
	Jd[1,1] = -l2*(Ad[0] + Ad[1])*np.sin(A[0] + A[1])
	return Jd

def compute_dynamics_terms(A, Ad, arm_params):
	a1,a2 = A[0], A[1]
	a1d,a2d = Ad[0], Ad[1]
	l1,l2 = arm_params['l'][0], arm_params['l'][1]
	m1,m2 = arm_params['m'][0], arm_params['m'][1]
	r1,r2 = arm_params['r'][0], arm_params['r'][1]
	i1,i2 = arm_params['i'][0], arm_params['i'][1]
	M = np.zeros((2,2))
	M[0,0] = i1 + i2 + (m1*r1*r1) + (m2*((l1*l1) + (r2*r2) + (2*l1*r2*np.cos(a2))))
	M[0,1] = i2 + (m2*((r2*r2) + (l1*r2*np.cos(a2))))
	M[1,0] = M[0,1]
	M[1,1] = i2 + (m2*r2*r2)
	C = np.zeros(2)
	C[0] = -(m2*l1*a2d*a2d*r2*np.sin(a2)) - (2*m2*l1*a1d*a2d*r2*np.sin(a2))
	C[1] = m2*l1*a1d*a1d*r2*np.sin(a2)
	return M,C

def inverse_dynamics(A,Ad,Add,arm_params):
	n = np.shape(A)[0]
	Q = np.zeros((n,2))
	for i in range(n):
	   M,C = compute_dynamics_terms(A[i,:],Ad[i,:],arm_params)
	   ACC = Add[i,:]
	   Q[i,:] = M.dot(ACC) + C
	return Q

def forward_dynamics(A0, Ad0, Q, t, arm_params):
	# Euler's method
	n   = np.shape(t)[0]
	A   = np.zeros((n,2))
	Ad  = np.zeros((n,2))
	Add = np.zeros((n,2))
	A[0,:], Ad[0,:], Add[0,:] = A0, Ad0, np.zeros(2)
	for i in range(n-1):
		M,C = compute_dynamics_terms(A[i,:], Ad[i,:], arm_params)
		J = arm_jacobian(A[i,:], arm_params)
		Add[i+1,:] = np.linalg.inv(M).dot(Q[i,:] - C)
		Ad[i+1,:] = Ad[i,:] + Add[i+1,:]*(t[i+1]-t[i])
		A[i+1,:] = A[i,:] + Ad[i+1,:]*(t[i+1]-t[i])
	return (A,Ad,Add)
