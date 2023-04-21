#include <traj.h>
#include <stdio.h>
#include <io.h>
#include <stdlib.h>
#include <ini.h>
#include <float.h>

int check_probabilities(double *probabilities, int prob_length){
    /**
    * routine that checks that input probabilities sum to 1
    *
    * Parameters
    * ----------
    *
    * `probabilities` : array of probabilities
    *
    * `prob_length` : array length
    */
   double sum = 0.0;
   int i;
   for (i = 0; i < prob_length; i++){
       sum += probabilities[i];
    }
   printf("sum of probs %lf \n", sum);
   if (sum > 0.99999 && sum < 1.00001){
       printf("properly normalized array of probabilities!\n");
       return 1;
    }
   else{
       printf("not normalized array...\n");
       return 0;
    }
}

void read_EnergyFile(char *EnergyFileName, traj *Trajectory){ 
    /**
    * routine that reads the input energy file
    *
    * Parameters
    * ----------
    *
    * `EnergyFileName` : energies filename
    *
    * `Trajectory` : traj object
    */
 
    FILE *fn;                   // energy file; 
    FILE *fe;                   // declaring error file 

    int rows_i, cols_i, i; 
    size_t line_buf_size = 0;

    char *token; //*string;  
    char *string = NULL; 

    /* Opening and checking if it exist and if it is empty */
    if(strcmp(EnergyFileName, "-1") == 0){
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. Energy file missing. It is required for the task chosen.\nAborting.\n");
        printf("Error. Energy file missing. It is required for the task chosen.\nAborting.\n");
        exit(EXIT_FAILURE);
    }
    else{
        fn = fopen(EnergyFileName, "r");

        printf("\nReading Energy FILE...\n");
        check_empty_file(fn, EnergyFileName);

        rows_i = n_rows(fn);       
        fseek(fn, 0, SEEK_SET);

        /* Checking if the number of rows is what we expect: number of rows == number of frames */
        if(rows_i > Trajectory->frames){
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. Energy file %s contains a number of rows(%d) higher than expected (%d frames)\n Aborting.\n", \
                    EnergyFileName, rows_i, Trajectory->frames); 
            fclose(fe); 
            printf("Error. Energy file %s contains a number of rows(%d) higher than expected (%d frames)\n Aborting.\n",\
                    EnergyFileName, rows_i, Trajectory->frames);  
            exit(EXIT_FAILURE);
        }

        if(rows_i < Trajectory ->frames){
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. Energy file %s contains a number of rows(%d) lower than expected (%d frames)\n Aborting.\n", \
                    EnergyFileName, rows_i, Trajectory->frames);
            printf("Error. Energy file %s contains a number of rows(%d) lower than expected (%d frames)\n Aborting.\n", \
                    EnergyFileName, rows_i, Trajectory->frames);
            exit(EXIT_FAILURE);
        }

        /* Checking for not-empty rows, only one column per row, each row is only FLOAT. Assign each float-string to energy[i] array */
        //string = (char*) malloc (rows_i);                                // Allocate char string 
	
	char out_filename[100];
        sprintf(out_filename, "./RESULTS/Probabilities_read.csv");
        FILE* f_prob;
        f_prob = open_file_w(out_filename);


        for(i=0; i<rows_i; i++) {

            getline(&string, &line_buf_size, fn);                    // Reading entire line. 

            check_empty_rows(string);                       	     // Checking for empty lines in your file. 

            token = strtok(string, " \t\v\r\n");                     // Splitting the string-line in columns separated by spaces, tabs, or \n  or \v or \r

            cols_i = columns(token);   		                     // Counting the number of columns in each row of file.  

            if(cols_i == 1){
                //check_float_string(string, i, EnergyFileName);       // If #column in each row is 1, then check that the column is a float number. 	//(!)
                Trajectory->energies[i] = atof(string);              // Assigning each float-string to energy[i] 1D-array. 
                Trajectory->energies_cg[i] = Trajectory->energies[i];              	//(!)
		fprintf(f_prob, "%.*e\n", 15, Trajectory->energies[i]); 			//(!)
            }

            if(cols_i > 1){
                fe = fopen("error.dat", "w"); 
                fprintf(fe,"Error. Each row must contain ONLY ONE column. In your file '%s', the %d-th row, has %d columns\n", EnergyFileName,i+1,cols_i);
                printf("Error. Each row must contain ONLY ONE column. In your file '%s', the %d-th row, has %d columns\n", EnergyFileName, i+1, cols_i);
                fclose(fe);
                exit(EXIT_FAILURE);
            }
        }
        fclose(fn);                                                       // Closing energy file. 
        free(string);
        //for(i=0; i<rows_i; i++)
        //    printf("energy[%d] = %f\n", i, Trajectory->energies[i]);
    }
}


void read_TrajectoryFile(char *TrajFileName, traj *Trajectory){
    /**
    * routine that reads the input xyz coordinate file
    *
    * Parameters
    * ----------
    *
    * `TrajFileName` : trajectory filename
    *
    * `Trajectory` : traj object
    */
 
    FILE *ft;                   // trajectory file; 
    FILE *fe;                   // declaring error file;  

    int rows_i, cols_i, i, p, frame_idx, j, count; 
    size_t line_buf_size = 0;

    char *token;
    char *token_arr[100];  //*string, *line;
    char *string = NULL; 
    char *line;  

    /* Opening and check if it exist and if it is empty */
    ft = fopen(TrajFileName, "r");
    printf("\nReading Trajectory FILE...\n");

    check_empty_file(ft, TrajFileName);

    rows_i = n_rows(ft);                                // Number of rows in trajectory file "ft".
    fseek(ft, 0, SEEK_SET);

    //string = (char *) malloc (rows_i/Trajectory->frames);                  // Allocate char string;
   
    line   = (char *) malloc (200); 					     // Allocate char line with size (lenght) = 200; 
									     // We are sure that the lenght of each line is less than 200

    /* Initialize the 2D-array Trajectory->traj_coords[][] */
    for(i = 0; i < Trajectory->frames; i++){
        for(j = 0; j < Trajectory->n_at + 1; j++){
            Trajectory->traj_coords[i][j] = 0.0;					//(!)
            /*Trajectory->traj_coords[i][3 * j + 0] = 0.0;
            Trajectory->traj_coords[i][3 * j + 1] = 0.0;
            Trajectory->traj_coords[i][3 * j + 2] = 0.0;*/
        }
    }

    /* Checking for empty rows; 
     * Checking that there are ONLY rows with one column and three columns 
     * Checking that the rows with one column correspond to an integer number i.e. the number of atoms; 
     * Checking that the rows with three columns are float corresponding to x,y,z trajectory coordinate. */

    frame_idx = 0;                                      // Initialize frame index to  0 
    p  = 0;                                             // Initialize counter "p" to 0 
    float sum = 0;    

    for(i=0; i<rows_i; i++){
        if(p>=Trajectory->n_at)                       // p increases from 0 to 3*atomnum i.e. p=[0;3*230) i.e. p = [0;689] //(!)
            p = 0; 

        getline(&string, &line_buf_size, ft);           // Reading entire line;   
   
        strcpy(line, string);                           // Copying "string" in "line"; we need it in case of three columns. 

        if( i != (Trajectory->n_at + 2)*(frame_idx-1) + 1)  // Checking for empty rows except for the 2nd row of each frame representing the title
            check_empty_rows(string);                   // that could also be an empty string   

        token = strtok(string, " \t\v\r\n");            // Splitting the string-line in columns separated by spaces, tabs, or \n  or \v or \r

        cols_i = columns(token);                        // Counting the number of columns in each row of file.  

         
        if(i!=(Trajectory->n_at + 2)*(frame_idx-1) + 1){   // exclude the 2nd row of each frame 

            if(cols_i == 1){
            
                if(i != (Trajectory->n_at + 2)*frame_idx) {
                    printf("Error. The %dth frame contains a different number of rows. Each frame must have %d rows\n", frame_idx+1, Trajectory->n_at + 2);
                    exit(EXIT_FAILURE);
                } 
                else{
                    frame_idx++;
		    if (i > 0) { //Saving the sum of all spins value for the frame [frame_idx - 2]
                        Trajectory->traj_coords[frame_idx - 2][Trajectory->n_at] = sum;
                        sum = 0;
                    }
                    if(frame_idx>Trajectory->frames){
                        fe = fopen("error.dat", "w");
                        fprintf(fe, "Error. The number of trajectory frames is higher than the declared one in parameter file (%d).\nAborting\n", Trajectory->frames);
                    	fclose(fe);
                    	printf("Error. The number of trajectory frames is higher than the declared one in parameter file (%d).\nAborting\n", Trajectory->frames);
                    	exit(EXIT_FAILURE);
                    }
                    check_int_string(string, i, TrajFileName);             // Checking that each row is an INTEGER number. 
                    if(atoi(string) != Trajectory->n_at){
                        fe = fopen("error.dat","w");
                    	fprintf(fe,"Error. The number of atoms at %dth row has length (%d) different from atomnum(%d). Aborting\n",\
                                i+1,atoi(string),Trajectory->n_at);
                   	 printf("Error. The number of atoms at %dth row has length (%d) different from atomnum(%d). Aborting\n",\
                                     i+1,atoi(string),Trajectory->n_at);
                    	fclose(fe); 
                    	exit(EXIT_FAILURE);
                    }
                }
            }
    
            /*if(cols_i == 2){								//(!)
 
                if(i != (Trajectory->n_at + 2)*frame_idx + 1){
                    fe = fopen("error.dat", "w");
                    fprintf(fe,"Error. Each row must not contain 2 columns (except the title in the 2nd row of each frame that can contain N columns).\n\
                      	        ONLY 1 or 4 columns are allowed in Trajectory (%dth row has 2 columns)\n", i+1);

                    printf("Error. Each row must not contain 2 columns (except the title in the 2nd row of each frame that can contain N columns)\n\
                            ONLY 1 or 4 columns are allowed in Trajectory (%dth row has 2 columns.\n", i+1);
                    fclose(fe); 
                    exit(EXIT_FAILURE);
                }
            }*/
            
            if (cols_i == 2) {								//(!)
                if (i != (Trajectory->n_at + 2) * frame_idx) {
                    token_arr[0] = strtok(line, " \t\v\r\n");

                    count = 0;
                    while (token_arr[count]) {
                        count++;
                        token_arr[count] = strtok(NULL, " \t\v\r\n");
                    }

                    for (j = 1; j <= count - 1; j++) {
                        //check_float_string(token_arr[j], i, TrajFileName);                     // Checking that each row is an FLOAT number.		//(!)
                        Trajectory->traj_coords[frame_idx - 1][p] = atof(token_arr[j]);          // Assigning each float-string to traj_coords[i][j] 2D-array 
			sum += Trajectory->traj_coords[frame_idx - 1][p];
                        p++;
                    }

                }

                else {
                    printf("Error. Each row of the frame (except the 1st row containing the number of atoms and the 2nd row containing the title),\n\
                    must contain 1 column corresponding to a spin state. Check also if there is one extra row in %dth frame\n",//CHANGED
                        frame_idx);
                    exit(EXIT_FAILURE);
                }
            }
            

            /*if(cols_i == 3){								//(!)
 
                if(i != (Trajectory->n_at + 2)*frame_idx + 1){
                    fe = fopen("error.dat", "w");
                    fprintf(fe,"Error. Each row must not contain 3 columns (except the title in the 2nd row of each frame that can contain N columns). \
                                ONLY 1 or 4 columns are allowed in Trajectory (%dth row has 3 columns)\n", i+1);

                    printf("Error. Each row must not contain 3 columns (except the title in the 2nd row of each frame that can contain N columns) \
                            ONLY 1 or 4 columns are allowed in Trajectory (%dth row has 3 columns.\n", i+1);
                    fclose(fe); 
                    exit(EXIT_FAILURE);
                }
            }


            if(cols_i == 4){ 
                if(i != (Trajectory->n_at + 2)*frame_idx){
                    token_arr[0] = strtok(line, " \t\v\r\n");

                    count = 0;  
                    while(token_arr[count]){
                        count++; 
                        token_arr[count] = strtok(NULL, " \t\v\r\n");
                    }

                    for (j=1; j <= count-1 ; j++){
                        check_float_string(token_arr[j], i, TrajFileName);                     // Checking that each row is an FLOAT number.
                        Trajectory->traj_coords[frame_idx-1][p] = atof(token_arr[j]);          // Assigning each float-string to traj_coords[i][j] 2D-array 
                        p++;
                    }

                }

                else{ 
                    printf("Error. Each row of the frame (except the 1st row containing the number of atoms and the 2nd row containing the title),\n\
                    must contain 4 columns corresponding to at_name, x, y, and z coodinates. Check also if there is one extra row in %dth frame\n",
                            frame_idx); 
                    exit(EXIT_FAILURE);
                }
            }*/

            if(cols_i > 2 ){ 								//(!)
                if(i != (Trajectory->n_at + 2)*frame_idx){ 
                    fe = fopen("error.dat","w"); 
                    fprintf(fe,     "Error. The maximum number of columns allowed is 2. The %dth row of your file contains %d columns\n", i+1, cols_i);
                    printf("Error. The maximum number of columns allowed is 2. The %dth row of your file contains %d columns\n", i+1, cols_i);
                    fclose(fe); 
                    exit(EXIT_FAILURE);
                }
            }
        }

    }  // END FOR LOOP 
    
    free(string);
    //free(token);
    free(line);
    fclose(ft);// Close trajectory file. 
    // 1st check: frame_idx must be = frames
    if (frame_idx !=  Trajectory->frames) {
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. Frames completed: %d, while declared frames in parameter file are %d. Trajectory incomplete.\nAborting\n",frame_idx,\
                     Trajectory->frames);
        fclose(fe);
        printf("Error. Frames completed: %d, while declared frames in parameter file are %d. Trajectory incomplete.\nAborting\n",frame_idx,\
                Trajectory->frames);
        exit(EXIT_FAILURE);
    }

    // final check: frame_idx corresponds with what we expect but the number of row is not the same => It means that 
    //              the frame completed are frame_idx - 1 
    else{
        if(rows_i !=  (Trajectory->n_at + 2)*Trajectory->frames){
            fe = fopen("error.dat", "w");
            fprintf(fe,"Error. Total number of rows: %d, while expected number of rows are %d. Trajectory incomplete.\nAborting\n",rows_i,\
                        (Trajectory->n_at + 2)*Trajectory->frames);
            fclose(fe);
            printf("Error. Total number of rows: %d, while expected number of rows are %d. Trajectory incomplete.\nAborting\n", rows_i,\
                    (Trajectory->n_at + 2)*Trajectory->frames);
            exit(EXIT_FAILURE);
        }
    }
}
