#define MAXCHAR  100				/* maximale Laenge der Pfadnamen */


/* ===================== return index for type of atom =============== */
int get_atom_type_index(char residuetype[MAXCHAR]);


/* ===================== return index for type of aa ================= */
int get_aa_index(char residuetype[MAXCHAR]);


/* ===================== return charge of atom ================= */
float return_charge(char residuetype[MAXCHAR]);
