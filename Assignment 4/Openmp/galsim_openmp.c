#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

int PRINT_DEBUG = 1; // 0 is off, 1 is light, 2 is intricate

typedef struct { // Storing the data for the particles in a Structure of Arrays (SoA)

    double *x; // x-positions
    double *y; // y-positions

    double *mass; // mass

    double *vx; // x-velocity
    double *vy; // y-velocity

    double *brightness;

} ParticleData;

// Function to print all particles (help from ChatGPT)
void print_particles(ParticleData particle_info, int N) {
    printf("\nParticle Data (x, y, mass, vx, vy, brightness):\n");
    for (int i = 0; i < N; ++i) {
        printf("Particle %d: %6f %6f %.6f %.6f %.6f %.6f\n",
               i, particle_info.x[i], particle_info.y[i], particle_info.mass[i],
               particle_info.vx[i], particle_info.vy[i], particle_info.brightness[i]);
    }
}

void free_particles_pointers(ParticleData *particles){ // Function for freeing all pointers in the ParticleData struct
    free(particles->x);
    free(particles->y);
    free(particles->mass);
    free(particles->vx);
    free(particles->vy);
    free(particles->brightness);
}

int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("Wrong number of input variables. Should follow format:\n./galsim N filename nsteps delta_t graphics");
        exit(1);
    }
    // Converitng the input to variables
    const int N = atoi(argv[1]);          // The number of stars/particles
    const char *filename = argv[2];       // The filename
    const int n_steps = atoi(argv[3]);    // The number of timesteps
    const double delta_t = atof(argv[4]); // The time step
    const int graphics = atoi(argv[5]);   // On or off, 1 or 0
    int n_threads = atoi(argv[6]);  // Number of threads in paralellization

    if (n_threads < 1){ // Checks for legal input for number of threads
        printf("Illegal number of threads, expected > 0");
        exit(1);
    }

    omp_set_num_threads(n_threads); // Setting the number of threads that the user inputted

    #pragma omp parallel num_threads(n_threads) // Checking that the parallelization works
    {
        int num_threads = omp_get_num_threads();
        printf("Thread %d out of %d\n", omp_get_thread_num(), num_threads);
        if (num_threads != n_threads){
            printf("Error in parallelization. Doesn't create the correct number of threads.\nWanted: %d   Created: %d", n_threads, num_threads);
            exit(1);
        
        }
    }

    // Reading the file and setting the values into the SoA

    ParticleData particles; // Initializes a particles SoA

    // Allocating memory for all arrays in the SoA of size N:
    particles.x = malloc(N * sizeof(double));
    particles.y = malloc(N * sizeof(double));
    particles.mass = malloc(N * sizeof(double));
    particles.vx = malloc(N * sizeof(double));
    particles.vy = malloc(N * sizeof(double));
    particles.brightness = malloc(N * sizeof(double));

    if (!particles.x || !particles.y || !particles.mass ||
        !particles.vx || !particles.vy || !particles.brightness){
            printf("Memory allocation failed");
            free_particles_pointers(&particles);
            exit(1);
        }
    
    FILE *file = fopen(filename, "rb");

    if (file == NULL){ // Checks that file got read
        printf("Error reading file");
        exit(1);
    }
    // File reading was done with the help of Chat-GPT
    for (int i = 0; i < N; ++i) {
        if (fread(&particles.x[i], sizeof(double), 1, file) != 1 ||
            fread(&particles.y[i], sizeof(double), 1, file) != 1 ||
            fread(&particles.mass[i], sizeof(double), 1, file) != 1 ||
            fread(&particles.vx[i], sizeof(double), 1, file) != 1 ||
            fread(&particles.vy[i], sizeof(double), 1, file) != 1 ||
            fread(&particles.brightness[i], sizeof(double), 1, file) != 1){
            printf("Error reading file at particle %d\n", i);
            fclose(file);

            // Freeing the pointers if failed
            free_particles_pointers(&particles);
            exit(1);
            
        }
    }

    fclose(file);


    // Beginning of the simulation, show the planets before the iteration
    if (PRINT_DEBUG > 0){
        print_particles(particles, N);
    }

    double const G = 100.0 / N; // Defining the instance dependent constant G

    // #######################
    // -----TIME STEPPING-----
    // #######################

    if (PRINT_DEBUG > 0){
        printf("Time stepping starts\n");
    }

    int chunk_size = N/n_threads;

    // Start time measurement
    double t_start = omp_get_wtime();

    for (int t = 0; t < n_steps; t++){ // Iterate through all time steps

        #pragma omp parallel for schedule(dynamic, 1) // Parallelizing the i loop with dynamic scheduling to combat load balancing problem
        for (int i = 0; i < N; i++){ // i is current particle

            // Reset force & acceleration vectors to {0, 0}:

            // i:s mass is loop invariant
            double i_mass = particles.mass[i];
            double r_hat_x = 0.0;
            double r_hat_y = 0.0;
            double F_x_i, F_y_i;
            double ax, ay;
	        double i_v_tot_x = 0.0;
	        double i_v_tot_y = 0.0;
            
            for (int j = i + 1; j < N; j++){ // Iterate through the particles that haven't been evaluated          

                // Doing the force calculation between particle i and j

                double r_x = particles.x[i] - particles.x[j]; // r vector x component
                double r_y = particles.y[i] - particles.y[j]; // r vector y component

                double r2 = r_x * r_x + r_y * r_y;
                double r = sqrt(r2);

                double r_eps_reci = 1 / (r + 0.001);
                
                r_hat_x = r_x * r_eps_reci; // Faster to multiply with reciprocal
                r_hat_y = r_y * r_eps_reci;

                double denom_rec = r_eps_reci*r_eps_reci; // = 1/(r + e_0)^2

                double F_scalar = particles.mass[j] * denom_rec;
                
                // Now that we have the total force, we can calculate the x- and y components
                F_x_i = F_scalar * r_hat_x;
                F_y_i = F_scalar * r_hat_y;                

                double F_x = -G * i_mass * F_x_i;
                double F_y = -G * i_mass * F_y_i;

                // Update acceleration for planet i using delta_t * (F/m)
                double i_mass_inv = 1 / i_mass;
                double j_mass_inv = 1 / particles.mass[j];

                // Acceleration contribution for i-j pair
                ax = (F_x * i_mass_inv);
                ay = (F_y * i_mass_inv);

                // Update planet i:s velocity using delta_t * acceleration
                i_v_tot_x += delta_t * ax;
                i_v_tot_y += delta_t * ay;
                
                // Planet j:s force is opposite
                
                particles.vx[j] += delta_t * (-F_x * j_mass_inv);                
                particles.vy[j] += delta_t * (-F_y * j_mass_inv);
                
            }

	    particles.vx[i] += i_v_tot_x;
	    particles.vy[i] += i_v_tot_y;

        }

        #pragma omp parallel for schedule(static, 1) // Parallelizes the for loop with static scheduling (even load)
        for (int i = 0; i < N; i++){

            // Update planet i:s position using delta_t * velocity
            particles.x[i] += delta_t * particles.vx[i];
            particles.y[i] += delta_t * particles.vy[i];
        }
        


    }

    // Time stepping done!
    
    // End time measurement
    double t_end = omp_get_wtime();
    double time_spent = t_end - t_start;

    if (PRINT_DEBUG > 0){
        printf("\nTime stepping done! Positions after simulation:\n");
        print_particles(particles, N);
        printf("Execution time: %.6f seconds\n", time_spent);
    }


    // - Writing output to file! (Helped by chat GPT-4o)

    FILE *output_file = fopen("result.gal", "wb");
    
    if (output_file == NULL){
        printf("File opening failed");
        fclose(file);
        free_particles_pointers(&particles); // Freeing the pointers
        exit(1);
    }

    for (int i = 0; i < N; ++i) {
        // Write the data for particle i in the correct order (x, y, mass, vx, vy, brightness, ax, ay)
        if (fwrite(&particles.x[i], sizeof(double), 1, output_file) != 1 ||
            fwrite(&particles.y[i], sizeof(double), 1, output_file) != 1 ||
            fwrite(&particles.mass[i], sizeof(double), 1, output_file) != 1 ||
            fwrite(&particles.vx[i], sizeof(double), 1, output_file) != 1 ||
            fwrite(&particles.vy[i], sizeof(double), 1, output_file) != 1 ||
            fwrite(&particles.brightness[i], sizeof(double), 1, output_file) != 1) {
            printf("Error writing to file");
            fclose(output_file);
            exit(1);}
            }
    

    fclose(output_file);

    // Freeing the pointers
    free_particles_pointers(&particles);

    // Program success
    
    return 0;
}
