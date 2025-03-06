#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

int PRINT_DEBUG = 1; // 0 is off, 1 is light, 2 is intricate

typedef struct { // Storing the data for the particles in a Structure of Arrays (SoA)

    double *x; // x-positions
    double *y; // y-positions

    double *mass; // mass

    double *vx; // x-velocity
    double *vy; // y-velocity

    double *brightness;

} ParticleData;

typedef struct {
    ParticleData *particles;
    int start;
    int end;
    int N;
    double delta_t;
    double G;
    pthread_mutex_t *mutex;
} ThreadedData;

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

void *threaded_function(void* arg){

    // Importing all the data
    ThreadedData *data = (ThreadedData*) arg;
    ParticleData *particles = data->particles;  // The information of all particles
    int start = data->start;                    // What planet to start from
    int end = data->end;                        // What planet to end on
    int N = data->N;                            // How many planets
    double delta_t = data->delta_t;             // The time step
    double const G = data->G;                   // The forceconstant (import or calculate?)

    // Starting

    for (int i = start; i < end; i++){ // i is current particle

            // Reset force & acceleration vectors to {0, 0}:

            // i:s mass is loop invariant
            double i_mass = particles->mass[i];
            double r_hat_x = 0.0;
            double r_hat_y = 0.0;
            double F_x_i, F_y_i;
            double ax, ay;

            for (int j = i + 1; j < N; j++){ // Iterate through the particles that haven't been evaluated                

                // Doing the force calculation between particle i and j

                double r_x = particles->x[i] - particles->x[j]; // r vector x component
                double r_y = particles->y[i] - particles->y[j]; // r vector y component

                double r2 = r_x * r_x + r_y * r_y;
                double r = sqrt(r2);

                double r_eps_reci = 1 / (r + 0.001);
                
                r_hat_x = r_x * r_eps_reci; // Faster to multiply with reciprocal
                r_hat_y = r_y * r_eps_reci;

                double denom_rec = r_eps_reci*r_eps_reci; // = 1/(r + e_0)^2

                double F_scalar = particles->mass[j] * denom_rec;
                
                // Now that we have the total force, we can calculate the x- and y components
                F_x_i = F_scalar * r_hat_x;
                F_y_i = F_scalar * r_hat_y;                

                double F_x = -G * i_mass * F_x_i;
                double F_y = -G * i_mass * F_y_i;

                // Update acceleration for planet i using delta_t * (F/m)
                double i_mass_inv = 1 / i_mass;
                double j_mass_inv = 1 / particles->mass[j];

                // Acceleration contribution for i-j pair
                ax = (F_x * i_mass_inv);
                ay = (F_y * i_mass_inv);

                // Update planet i:s velocity using delta_t * acceleration
                pthread_mutex_lock(data->mutex);
                particles->vx[i] += delta_t * ax;
                particles->vy[i] += delta_t * ay;
                
                // Planet j:s force is opposite
                particles->vx[j] += delta_t * (-F_x * j_mass_inv);
                particles->vy[j] += delta_t * (-F_y * j_mass_inv);
                pthread_mutex_unlock(data->mutex);
                
            }

        }
    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("Wrong number of input variables. Should follow format:\n./galsim N filename nsteps delta_t graphics n_threads");
        exit(1);
    }
    // Using converitng the input to variables
    // argv[0] = ./...)
    const int N = atoi(argv[1]);          // The number of stars/particles
    const char *filename = argv[2];       // The filename
    const int n_steps = atoi(argv[3]);    // The number of timesteps
    const double delta_t = atof(argv[4]); // The time step
    const int graphics = atoi(argv[5]);   // On or off, 1 or 0
    const int n_threads = atoi(argv[6]);  // The number of threads

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

    // Start time measurement
    clock_t start = clock();

    // Create the threads
    pthread_t threads[n_threads];
    ThreadedData threaded_data[n_threads];

    // Locking the variables?????
    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);

    //Calculating the total number of computations per thread
    int total_computation = (N * (N - 1)) / 2;
    double max_computation_per_thread = total_computation / n_threads;

    // Calculating where to start and end
    int start_comp[n_threads];
    int end_comp[n_threads];
    int count=0;
    start_comp[0] = 0;

    for (int j = 0; j < n_threads - 1; j++) {
        int sum = 0; // The sum is back to zero
        start_comp[j] = count;

        for (int i = count; i < N; i++) {
            sum += (N - i - 1); // Number of computations for particle i

            // Checking if the sum is bigger than the maximum sum
            if (sum >= max_computation_per_thread) {
                end_comp[j] = i;
                count = i + 1; // Next thread starts here
                break;
            }
        }
    }

    // We are at the end
    start_comp[n_threads-1] = end_comp[n_threads-2] + 1;
    end_comp[n_threads-1] = N - start_comp[n_threads-1];
    


    for (int t = 0; t < n_steps; t++){ // Iterate through all time steps

        for (int i=0; i<n_threads; i++){

            // Plugging the data into the thread
            threaded_data[i] = (ThreadedData) {
                .particles = &particles, 
                .start = start_comp[i],
                .end = end_comp[i]+1, // Ensures that the range includes enerything
                .N = N,
                .delta_t=delta_t,
                .G = G,
                .mutex = &mutex
                };
            pthread_create(&threads[i], NULL, threaded_function, &threaded_data[i]);
            
        }

        for (int i=0; i<n_threads; i++){

            // Joining all the threads
            pthread_join(threads[i], NULL);

        }


        for (int i = 0; i < N; i++){

            // Update planet i:s position using delta_t * velocity
            particles.x[i] += delta_t * particles.vx[i];
            particles.y[i] += delta_t * particles.vy[i];

        
            if (PRINT_DEBUG >= 2){
                printf("particle %d :s velocity is now %lf, %lf\n", i, particles.vx[i], particles.vy[i]);
            }
        }
        


    }

    // Time stepping done!
    
    // End time measurement
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

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
