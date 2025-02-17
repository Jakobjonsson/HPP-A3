#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int PRINT_DEBUG = 1;

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
        printf("Particle %d: %.3f %.3f %.3f %.3f %.3f %.3f\n",
               i, particle_info.x[i], particle_info.y[i], particle_info.mass[i],
               particle_info.vx[i], particle_info.vy[i], particle_info.brightness[i]);
    }
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Wrong number of input variables. Should follow format:\n./galsim N filename nsteps delta_t graphics");
        exit(1);
    }
    // Using converitng the input to variables
    // argv[0] = ./...)
    const int N = atoi(argv[1]);          // The number of stars/particles
    const char *filename = argv[2];       // The filename
    const int n_steps = atoi(argv[3]);    // The number of timesteps
    const double delta_t = atof(argv[4]); // The time step
    const int graphics = atoi(argv[5]);   // On or off, 1 or 0

    // Reading the file and setting the values into the SoA

    ParticleData particles; // Initializes a particles SoA

    FILE *file = fopen(filename, "rb");

    if (file == NULL){ // Checks that file got read
        printf("Error reading file");
        exit(1);
    }

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
            free(particles.x);
            free(particles.y);
            free(particles.mass);
            free(particles.vx);
            free(particles.vy);
            free(particles.brightness);
            exit(1);
        }

    for (int i = 0; i < N; i++) { // (Helped by chat GPT-4o)
        // Read the data for particle i in the correct order (x, y, mass, vx, vy, brightness)
        fread(&particles.x[i], sizeof(double), 1, file); // Position x
        fread(&particles.y[i], sizeof(double), 1, file); // Position y
        fread(&particles.mass[i], sizeof(double), 1, file); // mass
        fread(&particles.vx[i], sizeof(double), 1, file); // x velocity
        fread(&particles.vy[i], sizeof(double), 1, file); // y velocity
        fread(&particles.brightness[i], sizeof(double), 1, file); // Brightness
    }
    fclose(file);

    // Beginning of the simulation, show the planets before the iteration
    if (PRINT_DEBUG == 1){
        print_particles(particles, N);
    }

    double const G = 100/N; // Defining the instance dependent constant G

    /*
    -----TIME STEPPING-----
    */

    if (PRINT_DEBUG == 1){
        printf("Time stepping starts\n");
    }

    for (int t = 0; t < n_steps; t++){ // Iterate through all time steps

        for (int i = 0; i < N; i++){ // i is current particle

            // Reset force vector to {0, 0}:
            double F_x = 0.0, F_y = 0.0;

            // i:s mass is loop invariant
            double i_mass = particles.mass[i];

            for (int j = 0; j < N; j++){ // Iterate through all other particles to calculate total force exerted by the other particles

                if (i != j){ // Unless the current particle _is_ the other particle

                    // Doing the force calculation between particle i and j

                    double r_x = particles.x[i] - particles.x[j]; // r vector x component
                    double r_y = particles.y[i] - particles.y[j]; // r vector y component

                    double F_scalar_x = particles.mass[j] * ((r_x+0.001)*(r_x+0.001));
                    double F_scalar_y = particles.mass[j] * ((r_y+0.001)*(r_y+0.001));

                    // Now that the force between i and j is calculated, update the total force exerted on planet i

                    F_x += F_scalar_x;
                    F_y += F_scalar_y;

                }

            }
            F_x = -G * i_mass * F_x; // Doesn't need to be in loop (distributive)
            F_y = -G * i_mass * F_y; // Doesn't need to be in loop (distributive)

            // Calculate acceleration for planet i using the force divided by mass
            double i_mass_inv = 1 / i_mass;
            double ax = F_x * i_mass_inv;
            double ay = F_y * i_mass_inv;

            // Update planet i:s velocity using delta_t * acceleration
            particles.vx[i] += delta_t * ax;
            particles.vy[i] += delta_t * ay;

            if (PRINT_DEBUG == 1){
                printf("particle %d :s velocity is now %lf, %lf\n", i, particles.vx[i], particles.vy[i]);
            }

        }

        for (int i = 0; i < N; i++){

            // Update planet i:s position using delta_t * velocity
            particles.x[i] += delta_t * particles.vx[i];
            particles.y[i] += delta_t * particles.vy[i];
        }

    }

    // Time stepping done!

    if (PRINT_DEBUG == 1){
        printf("\nTime stepping done! Positions after simulation:\n");
        print_particles(particles, N);
    }


    // - Writing output to file! (Helped by chat GPT-4o)

    FILE *output_file = fopen("result.gal", "wb");
    for (int i = 0; i < N; ++i) {
        // Write the data for particle i in the correct order (x, y, mass, vx, vy, brightness)
        fwrite(&particles.x[i], sizeof(double), 1, output_file);
        fwrite(&particles.y[i], sizeof(double), 1, output_file);
        fwrite(&particles.mass[i], sizeof(double), 1, output_file);
        fwrite(&particles.vx[i], sizeof(double), 1, output_file);
        fwrite(&particles.vy[i], sizeof(double), 1, output_file);
        fwrite(&particles.brightness[i], sizeof(double), 1, output_file);
    }

    fclose(output_file);

    // Freeing the pointers
    free(particles.x);
    free(particles.y);
    free(particles.mass);
    free(particles.vx);
    free(particles.vy);
    free(particles.brightness);

    // Program success
    return 0;
}