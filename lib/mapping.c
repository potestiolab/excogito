/**
* \class cg_mapping_lib
* \brief library of functions that perform simple operations on CG mappings 
*/
//ciao
#include <stdio.h>
#include <stdlib.h>
#include <mapping.h>
#include <io.h>
#include <my_malloc.h>

void free_mapping(cg_mapping *mapping){
    /**
    * routine that frees the mapping
    * Parameters
    * ----------
    *
    * `mapping` : cg_mapping object
    */
    free_i1t(mapping->mapping);
    free_i1t(mapping->clusters);
    free_i1t(mapping->size);
    free_i1t(mapping->idx_cluster);//ADDED						//(!)
    free_i1t(mapping->omega);//ADDED							//(!)
    free(mapping);
}

void convert_mapping(cg_mapping *mapping, FILE *f_out) {
    /**
    * routine that prints out the mapping
    * Parameters
    * ----------
    *
    * `mapping` : cg_mapping object
    * 
    * `f_out` : file to write on
    */

    int i;
    fprintf(f_out, "conv mapping\n");
    for (i = 0; i < mapping->n_at; i++) {
        if (mapping->mapping[i] == 1) fprintf(f_out, "%d ", i);
    }
    fprintf(f_out, "\n");
}

void generate_random_mapping(cg_mapping *mapping, FILE *f_out) {
    /**
    * routine that generates a random mapping
    * 
    * Parameters
    * ----------
    *
    * `mapping` : cg_mapping object
    * 
    * `f_out` : file to write on
    */
 
    fprintf(f_out, "generating random mapping...\n");
    int i, j, rd;
    double x;
    for (i = 0; i < mapping->n_at; i++) mapping->mapping[i] = 0;
    j = 0;
    while (j < mapping->n_cg) {
        rd = rand() % mapping->n_at;
        if (mapping->mapping[rd] == 0) {
            mapping->mapping[rd] = 1;
            j++;
        }
    }
    fprintf(f_out, "random mapping generated!\n[");
    for (i = 0; i < mapping->n_at; i++) {
        fprintf(f_out, "%d ", mapping->mapping[i]);
    }
    fprintf(f_out, "]\n");
    // call convert mapping
    convert_mapping(mapping, f_out);
}

void update_mapping(cg_mapping *curr_mapping, cg_mapping *old_mapping, int frames){
    /**
    * routine that updates old_mapping with the data contained in curr_mapping
    * 
    * Parameters
    * ----------
    *
    * `curr_mapping` : current cg_mapping object
    * 
    * `old_mapping` : cg_mapping object to be updated
    * 
    * `frames` : length of the MD trajectory
    */
 
    int at_idx;
    // copying mapping
    for (at_idx = 0; at_idx < curr_mapping->n_at; at_idx++) {
        old_mapping->mapping[at_idx] = curr_mapping->mapping[at_idx];
    }
    // copying value of smap
    old_mapping->smap = curr_mapping->smap;
    // copying clusters and size
    int fr_idx;
    for (fr_idx = 0; fr_idx < frames; fr_idx++) {
        old_mapping->clusters[fr_idx] = curr_mapping->clusters[fr_idx];
        old_mapping->size[fr_idx] = curr_mapping->size[fr_idx];
        old_mapping->idx_cluster[fr_idx] = curr_mapping->idx_cluster[fr_idx];//ADDED	//(!)
        old_mapping->omega[fr_idx] = curr_mapping->omega[fr_idx];//ADDED		//(!)
    }
}

void read_MappingFile(char *MappingFileName, FILE *f_out_l, cg_mapping *mapping){
    /**
    * routine that reads the input mapping file
    *
    * Parameters
    * ----------
    *
    * `MappingFileName` : mapping filename
    *
    * `f_out_l` : output filename
    *
    * `cg_mapping` : cg_mapping object
    */
 
    FILE *fm;                   // mapping file; 
    FILE *fe;                   // declaring error file 

    int rows_i, cols_i, i, j;
    size_t line_buf_size = 0;

    int* mapp;                  // array mapping;

    char *token; //
    char *string = NULL; 

    /* Opening and checking if it exist and if it is empty */
    if(strcmp(MappingFileName, "-1") == 0){ 
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. Mapping file missing. It is required for the task chosen.\nAborting.\n");
        printf("Error. Mapping file missing. It is required for the task chosen.\nAborting.\n"); 
        exit(EXIT_FAILURE);
    }
    else{
        fm = fopen(MappingFileName, "r");
        printf("\nReading Mapping FILE...%s\n", MappingFileName);

        check_empty_file(fm, MappingFileName);

        rows_i = n_rows(fm);                            // Check for the number of rows in mapping file "fm".
        fseek(fm, 0, SEEK_SET);

        /* Checking if the number of rows is what we expect: number of rows == number of CG beads */
        if(rows_i > mapping->n_cg){
            fe = fopen("error.dat", "w");
            fprintf(fe,     "Error. Mapping file longer (%d atoms) than expected (%d). \n Aborting\n", rows_i, mapping->n_cg);
            printf("Error. Mapping file longer (%d atoms) than expected (%d). \n Aborting\n", rows_i, mapping->n_cg);
            exit(EXIT_FAILURE); 
        }

        if(rows_i < mapping->n_cg){
            fe = fopen("error.dat", "w"); 
            fprintf(fe,     "Error. Mapping file shorter (%d atoms) than expected (%d). \n Aborting\n", rows_i, mapping->n_cg);
            printf("Error. Mapping file shorter (%d atoms) than expected (%d). \n Aborting\n", rows_i, mapping->n_cg);
            exit(EXIT_FAILURE);
        }

        /* Checking for empty row; Checking for only one column per row; Checking that each row is an INTEGER number;  
         * Checking that each row is a number between 0 and the number of atoms in atomistic structure, i.e. [0, n_at) */

        // string = (char*)malloc(rows_i);                                // Allocate char string. 
        mapp   = (int*)malloc(sizeof(int)*rows_i);                     // Allocate mapp[i] 1D-array. 

        for(i=0; i<rows_i; i++){

            getline(&string, &line_buf_size, fm);      		       // Reading entire line; 

            check_empty_rows(string);                   	       // Checking for empty rows 

            token = strtok(string, " \t\v\r\n");        	       // Splitting the string-line in columns separated by spaces, tabs, or \n  or \v or \r 

            cols_i = columns(token);                                   // Counting the number of columns in each row of file.  

            if(cols_i == 1){
                check_int_string(string, i, MappingFileName);          // Checking that each row is an INTEGER number. 
                mapp[i] = atoi(string);                                // Assigning each integer-string to mapp[i] 1D-array

                if(mapp[i] < 0 || mapp[i] >= mapping->n_at){           // Checking that each row has an index between 0 and N_atoms (not included) 
                    fe = fopen("error.dat", "w");
                    fprintf(fe, "Error. Mapping file contains atom ID %d at %dth row, while declared number of atoms is %d. The range is [0, n_atoms)\n \
                                 Aborting.\n", mapp[i], rows_i+1, mapping->n_at);
                    printf("Error. Mapping file contains atom ID %d at %dth row, while declared number of atoms is %d. The range is [0, n_atoms)\n \
                                 Aborting.\n", mapp[i], rows_i+1, mapping->n_at);
                    fclose(fe); 
                    exit(EXIT_FAILURE);
                }
            }

            if(cols_i > 1){
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. Mapping file must contain ONLY ONE column per row. In your file '%s' the %dth row, has %d columns\n", \
                             MappingFileName, rows_i+1, cols_i);
                printf("Error. Mapping file must contain ONLY ONE column per row. In your file '%s' the %dth row, has %d columns\n", \
                                 MappingFileName, rows_i+1, cols_i);
                fclose(fe);
                exit(EXIT_FAILURE);
            }
        }

        /* Checking for duplicates */
        for(i = 0; i < mapping->n_cg ; i++) {
            for(j = i + 1; j < mapping->n_cg; j++) {
                if(mapp[i] == mapp[j]){
                    fe = fopen("error.dat", "w");
                    fprintf(fe, "Error. Mapping file contains duplicate atom ID %d.\n Aborting\n", mapp[i]); 
                    printf("Error. Mapping file contains duplicate atom ID %d.\n Aborting\n", mapp[i]);
                    fclose(fe); 
                    exit(EXIT_FAILURE);
                }
            }
        }

        fclose(fm);
        } // END first ELSE in fm FILE

        /* Filling out binary array mapping->mapping[i]. It has value 1 if AT_index is found, otherwise it is Zero. */
        for(i=0; i<mapping->n_cg; i++){
            for(j=0; j<mapping->n_at; j++){
                if(mapp[i] == j)
                    mapping->mapping[j] = 1; 
            }
        }

        for(i=0; i<rows_i; i++)
            printf("mapping[%d] = %d\n", i, mapp[i]);

        fprintf(f_out_l,"loaded mapping\n[");

        for (i = 0; i < mapping->n_at; i++) 
            fprintf(f_out_l,"%d ", mapping->mapping[i]);
        
        fprintf(f_out_l,"]\n");
        convert_mapping(mapping,f_out_l);

        free(string); 
        free(mapp); 
}

/*void load_mapping_matrix(char *mappings_filename, FILE *f_out_l, cg_mapping *mapping_matrix[], int nmaps) {*/
    /** 
    * routine that reads the input mapping matrix 
    * 
    * Parameters
    * ----------
    *
    * filename : mapping filename
    * 
    * f_out_l : output filename
    * 
    * cg_mapping : cg_mapping object
    */
   /*
   printf("loading %d mappings\n", nmaps);
   printf("with cgnum = %d\n",mapping_matrix[0]->n_cg);
   char delim[] = " ";
   FILE* filePointer;
   int bufferLength = 5000;
   char buffer[bufferLength];
   filePointer = fopen(mappings_filename, "r");
   int map_idx = 0; // index of mapping
   int at_idx; // atomic index
   int knt; // counter over spaces
   while(fgets(buffer, bufferLength, filePointer)) {
        knt = 0;
        printf("%s\n", buffer);
        char *ptr = strtok(buffer, delim);
        while(ptr != NULL){
            if (knt < mapping_matrix[map_idx]->n_cg){
                at_idx = atoi(ptr);
                mapping_matrix[map_idx]->mapping[at_idx] = 1;
            }
            ptr = strtok(NULL, delim);
            knt += 1;
        }
        map_idx +=1;
   }
   fclose(filePointer);
}*/

void read_mapping_matrix(char *mappings_filename, FILE *f_out_l, cg_mapping *mapping_matrix[], int nmaps) {
    /** 
    * routine that reads the input mapping matrix 
    * 
    * Parameters
    * ----------
    *
    * `filename` : mapping filename
    * 
    * `f_out_l` : output filename
    * 
    * `cg_mapping` : cg_mapping object
    * 
    * `nmaps` : number of mappings defined in parameter file
    */
 
    FILE *fm;                   // Mapping-Matrix FILE; 
    FILE *fe;                   // Declaring error FILE;  

    int rows_i, cols_i, i, j, k, p;
    size_t line_buf_size = 0;

    int map_idx = 0;            // index of mapping
    int at_idx;		       	// atomic index

    int **matrix;               // 2D array mapping matrix 

    char *token, *string, *line;

    printf("loading %d mappings\n", nmaps);
    printf("with cgnum = %d\n", mapping_matrix[0]->n_cg);

    /* Opening and checking if it exist and if it is empty */
    fm = fopen(mappings_filename, "r");
    printf("\nReading Mapping FILE...%s\n", mappings_filename);

    check_empty_file(fm, mappings_filename);

    rows_i = n_rows(fm);                                                 // Check for the number of rows in mapping file "fm".
    fseek(fm, 0, SEEK_SET);

    matrix = malloc(nmaps*sizeof(int*));                                 // Allocate 2D array matrix[nmps][n_cg]
    for (i=0; i<nmaps; i++){
        matrix[i] = malloc(mapping_matrix[0]->n_cg*sizeof(int)); 
    }   

    /* Checking if the number of rows is what we expect: number of rows == number of mappings (n_mappings in parameter file) */
    if(rows_i > nmaps){
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. Mapping Matrix file longer (%d rows) than expected (%d).\nAborting\n", rows_i, nmaps);
        printf("Error. Mapping file longer (%d rows) than expected (%d).\nAborting\n", rows_i, nmaps);
        exit(EXIT_FAILURE); 
    }

    if(rows_i < nmaps){
        fe = fopen("error.dat", "w"); 
        fprintf(fe,"Error. Mapping Matrix file shorter (%d rows) than expected (%d).\nAborting\n",rows_i, nmaps);
        printf("Error. Mapping Matrix file shorter (%d rows) than expected (%d).\nAborting\n", rows_i, nmaps);
        exit(EXIT_FAILURE);
    }

    /* Checking also that the number of rows (n_mappings) is shorter than 5000, otherwise the output could be too big */ 
    if(rows_i > 5000){
        fe = fopen("error.dat", "w"); 
        fprintf(fe, "Error. Mapping Matrix file is too big. Please use a file shorter(%d rows) than 5000 rows\n", rows_i); 
        printf("Error. Mapping Matrix file is too big. Please use a file shorter(%d rows) than 5000 rows\n", rows_i);
        exit(EXIT_FAILURE);  
    }

    /* Checking for empty rows; Checking that each row contains a number of columns == n_cg columns;
     * Checking that each element (at_idx) is an INTEGER number;  
     * Checking that each element is a number between 0 and the number of atoms in atomistic structure, i.e. [0, n_at) */

    string = (char*)malloc(5*mapping_matrix[0]->n_cg);                                	// Allocate char string. 
    line   = (char*)malloc(5*mapping_matrix[0]->n_cg);                                	// Allocate char line; 

    for(map_idx=0; map_idx<rows_i; map_idx++){

        getline(&string, &line_buf_size, fm);      		   	// Reading entire line; 

        strcpy(line, string);                                      	// Copying "string" in "line"; we need it in case of "n_cg" columns. 

        check_empty_rows(string); 	                  	        // Checking for empty rows 

        token = strtok(string, " \t\v\r\n");        	                // Splitting the string-line in columns separated by spaces, tabs, or \n  or \v or \r 
        cols_i = columns(token);                                        // Counting the number of columns in each row of file.  
    
        if(cols_i == mapping_matrix[0]->n_cg){
       
            token = strtok(line, " \t\v\r\n");          
            p = 0; 
            while(token){
                check_int_string(token, map_idx, mappings_filename);       // Checking that each element is an INTEGER number.
                at_idx = atoi(token);
                matrix[map_idx][p] = atoi(token);                          // matrix[][]: it needs for checking duplicates.  
                p++; 
  
                if(at_idx< 0 || at_idx >= mapping_matrix[0]->n_at){        // Checking that each element has an index between 0 and N_atoms (not included) 
                    fe = fopen("error.dat", "w");                          // \r brings your cursor to the beginning of the line
                    fprintf(fe,"Error. Mapping Matrix file contains atom ID %d at %dth row, while declared number of atoms is %d.\
                               \rThe range is [0, n_atoms).\nAborting.\n", at_idx, map_idx+1, mapping_matrix[0]->n_at);
                    printf("Error. Mapping Matrix file contains atom ID %d at %dth row, while declared number of atoms is %d.\
                           \rThe range is [0, n_atoms).\nAborting.\n", at_idx, map_idx+1, mapping_matrix[0]->n_at);
                    fclose(fe); 
                    exit(EXIT_FAILURE);
                }
                
                mapping_matrix[map_idx]->mapping[at_idx] = 1;
                token = strtok(NULL, " \t\v\r\n");
            }
        } 

        if(cols_i != mapping_matrix[0]->n_cg){
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. Each line of Mapping Matrix file must contain %d (n_cg) columns per row.\n\
                        \rIn your file '%s' the %dth row, has %d columns.\n", mapping_matrix[0]->n_cg, mappings_filename, map_idx+1, cols_i);
            printf("Error. Each line of Mapping Matrix file must contain %d (n_cg) columns per row.\n\
                   \rIn your file '%s' the %dth row, has %d columns.\n", mapping_matrix[0]->n_cg, mappings_filename, map_idx+1, cols_i);
            fclose(fe);
            exit(EXIT_FAILURE);
        }
       
        /* Checking for duplicates */ 
        for(i = 0; i < mapping_matrix[0]->n_cg; i++){
            for(j = i+1; j < mapping_matrix[0]->n_cg; j++){

                if(matrix[map_idx][i] == matrix[map_idx][j]){
                    fe = fopen("error.dat", "w");
                    fprintf(fe,"Error. Mapping file contains duplicate atom ID %d in the mapping of %dth row.\n Aborting\n",matrix[map_idx][i],map_idx+1);
                    printf("Error. Mapping file contains duplicate atom ID %d in the mapping of %dth row.\n Aborting\n", matrix[map_idx][i],map_idx+1);
                    fclose(fe); 
                    exit(EXIT_FAILURE);
                 }
            }
        }
    
    }   //END first FOR
    fclose(fm);
} 
