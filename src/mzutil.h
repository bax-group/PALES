#define	MAXCHAR  100
#define MAX_OPTIONS 10
#define Real  double


/* standard error handler */
void nerror(char error_text[]);

/* ------- error handler for missing file ------- */
void nerror_datei(char error_txt[]);

/* ------- Kontrollmeldung fuer Einlesen einer Datei ------- */
void KontrollMeldung(char error_txt[]);

/* ------- Kontrollmeldung fuer Einlesen einer Datei ------- */
void KontrollMeldung_Lesen(char error_txt[]);

/* ------- Kontrolle, ob Programm richtig gestartet wurde ------- */
void Aufruf_Check(int noS_files, int noI_files, void (*error_out)(char *));

/* ------- dot product ------- */
Real dot_product(Real *vector1, Real *vector2);

/* ------- norm of vector ------- */
Real norm_vector(Real *vector);


int Gleich(char *L,char *Mit);

void Matrix_Initialisierung(Real **Matrix, int row_offset, int row_length,
			 int col_offset, int col_length, Real initial_wert);

/* ------- print matrix ------- */
void print_matrix(Real **Matrix, int no_row, int no_col);
