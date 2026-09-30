#define MAX_M 400			/*  # of rows, max # of dipolar couplings  */
#define MAX_NO_DC 2000
#define	MAX_ATOM_NO 5000
#define	MAX_NO_PT 5000
#define	MAX_NO_PT_SQR 5000
#define	MAX_ORIENT_PT 1000
#define	RADIUS 1.0
#define Real  double
#define MAXCHAR  100				/* maximale Laenge der Pfadnamen */
#define DC_THRESH  1000				/* cutoff for allowed DCs */
#define MAX_DIG_POT 10000
#define MAX_DOTS 3000
#define MAX_ENSEMBLE 8000


#define	MAXCHAR  100
#define MAX_OPTIONS 10

#define FABS(X) (((X) < 0.0) ? -(X) : (X))

#define NO_LIST 4
#define H_NO 329
#define N_ind 211
#define HN_ind 551
#define CA_ind 114
#define HA_ind 476
#define C_ind 1


Real nherror=1.5;


int debug=0;			/* debugging switched off */


int  liquid_type=1;		/* bicelle */
Real bic_phage=0;		/* take rim of bicelle into account */
int  grid_type=2;		/* 1D grid */
Real concentration=0.05;	/* total lipid concentration [%] */
Real ellipticity=1.0;		/* not used */
Real max_obs_dim=20.0;		/* radius of bicelle rim or of phage (30.0) [A] */
Real grid_spacing=0.5;		/* spacing of grid [A] */
Real medium_dist_max=8000.0;	/* not used */
Real molar_ratio=3.0;		/* DMPC/DHPC molar ratio */
Real order_parameter=0.8;	/* LIQUID_CRYSTALL_ORDER: [0-1] */
int  select_flag=0;		/* complete protein is taken into account (1/2/12/3) */
int  refine_ok=1,scale_down_ok=1; /* REFINE_SAUPE: & SCALE_SAUPE: */
int  mode=0;			/* newest program mode (0/1/2)*/
int  H_flag=0;			/* flag for taking protons into account (0/1) */
Real rotate_ok=0;		/* flag for allowing rotation of protein */
Real atom_radius=3.0;		/* radius of pseudo atoms [A] */
Real radius_bic=200, order_rel=0.25;
Real step_out_length;


/*** output ***/
int pr_eigen=0;			/* flag for printing eigensystem (0/1) */
int svd_flag=1;			/* flag for doing SVD */
int coup_flag=1;
int charge_info=1;


/*** electrostatics ***/
int electrostatics=0;		/* flag for incorporation of electrostatics (0/1)*/
int pKa_default=0;		/* using pKa default values */
int ch_termini=1;		/* flag for putting charges at termini (0/1)*/
double k_ener=1.0;		/* scaling for EL energy */
double k_pot=1.0;		/* scaling for EL potential */
float ch_weight=1.0;		/* weight charges */
float ch_neg_scale=1.0;		/* scaling down of negative charges */
float ch_pos_scale=1.0;		/* scaling down of positive charges */
float pH_auto=5.4;
float cutoffPOT = 100.0;	/* cutoff for taking into account */
				/* any interaction (Angstrom) */
float salt_conc = 0.05;		/* [M] molarity mol/l */
float surf_ch_dens = 0.8;	/* surface charge density [e/nm2] */
float temperaturePot =  298.150;	/* temperature [K] */


/*** digitization of spherical grid ***/
int old_mode = 0;
int dot_density = 100;
int dig_theta=9,dig_phi=18,dig_psi=18;


/*** ensemble averaging ***/
int	ensemble_out=0;		/* put out surface ensemble */
int	ensemble_averaging=0;	/* ensemle averaging flag (0/1) */
int	ensemble_density=200;	/* density of dots on ensemble sphere */
int     change_angle=0;		/* change rel. orientation of domains (0/1) */
float   alphaDom2=0,betaDom2=0,gammaDom2=0;  /* Euler angles for adjusting */
					     /* angle between 2 domains */
int     resDom1Origin=0; 	/* origin for orientation of domain 1 */
int	resDom2Origin=0; 	/* origin for orientation of domain 2 */
int	startResDom2,startResDom1,endResDom1,endResDom2;
				/* domain boundaries */
int     digit_ensemble=1; 	/* digitization of rotation about */
				/* 3rd Euler angle */
/* limits for polar angles of cone */
float	betaConeL = 0.0;
float   betaConeH = 30.0;
float   alphaConeC = 0.0;
float   alphaConeErr = 180.0;
float   alphaConeCInv = 0.0;

int     startResPart1,endResPart1,startResPart2,endResPart2;




/* ========================================================================= */
/* ========================================================================= */
/* ========================================================================= */


/* ===================== calculate center of mass ========================= */
Real calculate_center_mass(int no_of_masses, Real *x, Real *y,
	Real *z, Real *mass, Real *center_of_mass);


/* ===================== shift coordinate system =========================== */
void shift_coordinate_system(Real *new_origin, int no_atoms, 
	Real *x, Real *y, Real *z);


/* ============== calculate maximum diameter of protein ==================== */
void calculate_extreme_length(int no_atoms, Real *x, Real *y, Real *z,
	Real *max_distance, Real *min_distance);


/* ===================== make grid =================================== */
int make_grid(Real s_l_ratio, Real min_protein_length, Real max_protein_length,
	Real max_obstruction, int grid_flag, int liquid_flag,
	Real spacing, Real grid[MAX_NO_PT_SQR][4], int *no_pt_ph,
	Real grid_ph[MAX_NO_PT_SQR][4], Real *final_distance,
	Real *final_distance_ph);


/* ===================== setup model =================================== */
int setup_model(Real conc_ratio, Real s_l_ratio, Real min_protein_length,
	Real max_protein_length, Real max_obstruction, int liquid_flag,
	int grid_flag, Real spacing, Real grid[MAX_NO_PT_SQR][4], 
	int *no_pt_ph, Real grid_ph[MAX_NO_PT_SQR][4],
	Real *final_distance, Real *final_distance_ph);


/* ======================== rotate vector ================================= */
void rotate_vector(Real theta, Real phi, Real psi, Real *input_vector,
	Real *rotated_vector);


/* ================== rotate coordinate system =========================== */
void rotate_coorsystem(Real theta, Real phi, Real psi, Real *x_axis_tensor,
	Real *y_axis_tensor, Real *z_axis_tensor);


/* ======= calculation of Szz according to Prestegard ===================== */
Real Szz(Real tensor_zaxis[4]);


/* ========== calculation of Sxx according to Prestegard =============== */
Real Sxx(Real tensor_xaxis[4]);


/* ========== calculation of Syy according to Prestegard =============== */
Real Syy(Real tensor_yaxis[4]);


/* ============= calculation of Sxy according to Prestegard ============== */
Real Sxy(Real tensor_xaxis[4], Real tensor_yaxis[4]);


/* ============ calculation of Sxz according to Prestegard =============== */
Real Sxz(Real tensor_xaxis[4], Real tensor_zaxis[4]);


/* ============ calculation of Syz according to Prestegard =============== */
Real Syz(Real tensor_yaxis[4], Real tensor_zaxis[4]);


/* ============== simulate sterically induced alignment =================== */
void ssia(Real max_obstruction, int no_atoms, Real *x, Real *y, Real *z,
	double max_protein_length, double min_protein_length,
	int digit_psi, int n_dots, double *dots, Real saupe[5]);


/* ================ simulate macromolecular alignment ====================== */
void sma(Real max_obstruction, int no_atoms, Real *x, Real *y, Real *z,
	int digit_psi, int n_dots, double *dots,
	double max_protein_length, double min_protein_length,
	int no_of_ch, int *chargeAtomID,
	float *charges, int digPotential, Real potential[MAX_DIG_POT],
	Real offset, Real stepsize, Real saupe[5]);


/* ===================== substract two matrices ========================= */
void print_alignment(Real psiList[8],Real thetaList[8],Real phiList[8],
	Real da, Real dr);

/* ============ select atoms which should be considered in simulation ============= */
int select_atoms(int tot_no_atoms, Real *x, Real *y, Real *z, 
	Real *mass, int *atom, int *aa, int *atomno, int *res_nr);


/* ============ select atoms which should be considered in simulation ============= */
void calculate_scaling(Real final_distance, Real final_distance_ph,
	Real max_obstruction, int *no_pt, int *no_pt_ph);


/* ===================== return mass for type of atom ================ */
float return_mass(char residuetype[MAXCHAR]);


/* ========== compare measured with predicted dipolar couplings ======== */
void comp_meas_pred(Real saupe[5], int nodc, int *residue, Real *dc,
	int *dcmax_ind, Real *coor_x, Real *coor_y, Real *coor_z);


/* ================== predict dipolar coupling ==================== */
Real findDmz( Real xIJ_temp, Real yIJ_temp, Real zIJ_temp, 
           Real diVal, Real saupe[5]); 


/* ===== calculation of DCs from Saupe matrix in principal order frame ==== */
Real dc_calc(int type, Real coor_sys[3][3], Real prot_vec[3],
	Real tensor[3]);


/* ============= take into account DHPC in solution ============ */
void refine_bicelle();


/* ========================================================================= */
int compareDia(Real eigVal[3], Real eigVec[3][3], int no_dcs, int *res,
		Real *dc, Real *coor_x, Real *coor_y, Real *coor_z,
		Real predicted[16][MAX_ATOM_NO], float *rmsdList, float *r_statList,
		float *prob_statList, float *z_statList);


/* ========================================================================= */
int predict_comp_dia(Real eigVal[3], Real eigVec[3][3], int no_dcs, int *res,
		Real *dc, Real *coor_x, Real *coor_y, Real *coor_z,
		Real predicted[16][MAX_ATOM_NO], float *rmsdList, float *r_statList,
		float *prob_statList, float *z_statList);


/* ========================================================================= */
int make_eigVecList(Real eigVec[3][3], Real eigVecXlist[16][3],
	Real eigVecYlist[16][3], Real eigVecZlist[16][3]);


/* =============================================================== */
int cross_productM2(int n, int m, Real eigVec[3][3], Real product[3]);


/* =============================================================== */
int add_eigVec_List(Real eigVecXlist[16][3], Real eigVecYlist[16][3],
	Real eigVecZlist[16][3], int n, Real eigVec[3][3]);


/* =============================================================== */
int mirror_solutions(int n, Real eigVec[3][3], Real eigVecXlist[16][3],
	Real eigVecYlist[16][3], Real eigVecZlist[16][3]);


/* =============================================================== */
void print_compareDia(int n, int no_dcs, int *res, Real *dc, 
		Real predicted[16][MAX_ATOM_NO], float *rmsdList,
		float *r_statList, float *prob_statList, float *z_statList);


/* =============================================================== */
int setup_svd(int no_dcs, Real *coor_x, Real *coor_y, Real *coor_z,
	Real *dc, int *dcmax_ind, float *x1, float *y1, float *z1,
	float *x2, float *y2, float *z2,
	float *dc2, float *dcmax, float *error);


/* ===============    removal of protons   ====================== */
int remove_protons(int n, Real *x, Real *y, Real *z, 
	Real *mass, int *atom, int *aa, int *atomno, int *res_nr);


/* ===============    remove residues   ====================== */
int remove_residues(int n, Real *x, Real *y, Real *z, 
	Real *mass, int *atom, int *aa, int *atomno, int *res_nr,
	int no_select, int select[500]);


/* ===============    select atom types  ====================== */
int select_atomtypes(int n, Real *x, Real *y, Real *z, 
	Real *mass, int *atom, int *aa, int *atomno, int *res_nr,
	int no_select, int select[500]);


/* ===============    select atom types  ====================== */
void  legal_notice();


/* =============== */
void print_eigensystem(int flag, Real saupe[5], Real eigVal[3], Real eigVec[3][3]);


/* ===============  determine deviation of axes  ====================== */
int calc_dcos(Real eigVec [3][3], Real eigVec_svd[3][3]);

/* ============   select atoms taking care of charges  ================= */
int do_atom_selection(int tot_no_atoms, Real *x, Real *y, Real *z, 
	Real *mass, int *atom_type, int *aa_type, int *atomnumber,
	int *residue_nr, int no_of_ch, int *chargeAtomID);

/* =============== */
/* =============== */
/* =============== */


/* ===================== open a file ============================= */
FILE *DateiOeffnen(char dateiname[MAXCHAR],char *access);

/* ===================== reading pdb file ============================= */
int reading_PDBfile(char dateiname[MAXCHAR], Real *x, Real *y, Real *z,
	Real *mass, int *atom, int *aa, int *atomno, int *res_nr);

/* ================ reading measured dipolar couplings =================== */
int reading_DCfile(int flags[MAX_OPTIONS], FILE *flag_file[MAX_OPTIONS],
	int *residue, Real *dc, int *dcmax, int *dab_sign);

/* ===================== extraction of vectors ============================= */
void reading_vectors(int nodc, int *residue, int *dcmax, int noatoms, Real *x,
	Real *y, Real *z, int *atom, int *res_nr,
	Real *coor_x, Real *coor_y, Real *coor_z);

/* =========================== reading pKa values ======================= */
int read_flags(int inputc, char *input[], int *flag_ok, int flags[MAX_OPTIONS],
	FILE *flag_file[MAX_OPTIONS]);

/* ===================== read configuration file ============================ */
void read_auto(int flags[MAX_OPTIONS], FILE *flag_file[MAX_OPTIONS]);

/* ===================== reading pKa file ============================ */
int reading_pKafile(int flags[MAX_OPTIONS], FILE *flag_file[MAX_OPTIONS],
	int noatoms, int *atom, int *aa, int *res_nr,
	int *chargeAtomID, float *charges);

/* ===================== reading potential ============================ */
int reading_Potential(int flags[MAX_OPTIONS],FILE *flag_file[MAX_OPTIONS],
	Real potential[MAX_DIG_POT], Real *offset, Real *stepsize);

/* ========================================================================= */
/* ========================================================================= */

/* ============================= rotate protein =========================== */
void  rotate_protein(double angle[3], int no_atoms, Real *x, Real *y,
	Real *z);

/* ====================== rotate part of molecule ========================== */
void  rotate_domain(double angle[3], int noatoms, Real *x, Real *y, Real *z,
	int *atom, int *aa, int *atomno, int *res_nr, 
	int start_res, int end_res);

/* ========== compare measured with predicted dipolar couplings ======== */
void statistics_data(int flags[MAX_OPTIONS],FILE *flag_file[MAX_OPTIONS],
	int nodc, int *residue, float *predicted, Real *dc,
	int *dcmax_ind);

/* ========================================================================= */
void init_orient(int noatoms, double *x, double *y, double *z, 
	Real *mass, int *atom, int *aa, int *atomno, int *res_nr);

/* ========================================================================= */
void postitionDomain(int noatoms, double *x, double *y, double *z, 
	int *atom, int *aa, int *atomno, int *res_nr,
	int origin1, int origin2, int start_res, int end_res);

/* ========================================================================= */
void  calculate_inertia(int noatoms, Real *x, Real *y, Real *z,
	int *atom, int *aa, int *atomno, int *res_nr, 
	int start_res, int end_res, double *inertiaTensor);

/* ========================================================================= */
int writeRot_PDBfile(char dateiname[MAXCHAR], int noatoms, 
	Real *x, Real *y, Real *z);


/* ========================================================================= */
void readEnsemble(int flags[MAX_OPTIONS],FILE *flag_file[MAX_OPTIONS],
	int *n, double *alpha, double *beta, double *gamma);

/* ========================================================================= */
void comp_meas_predEn(FILE *fp, Real saupe[5], int nodc, 
	int *residue, Real *dc, int *dcmax_ind, Real *coor_x, 
	Real *coor_y, Real *coor_z);

/* ========================================================================= */
int calc_dcosEn(Real eigVec [3][3], Real eigVec_svd[3][3], float dcos[3]);
