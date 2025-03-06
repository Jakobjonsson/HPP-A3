# include <stdio.h>
# include <stdlib.h>

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
        printf("Particle %d: %.6f %.6f %.6f %.6f %.6f %.6f\n",
               i, particle_info.x[i], particle_info.y[i], particle_info.mass[i],
               particle_info.vx[i], particle_info.vy[i], particle_info.brightness[i]);
    }
}

int main(int argc, char *argv[]){
    if (argc != 3){
        printf("Failed. Provide input path");
        exit(1);

    }
    const int N = atoi(argv[1]);
    const char *filename = argv[2]; // The filename

    
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
            free(particles.x);
            free(particles.y);
            free(particles.mass);
            free(particles.vx);
            free(particles.vy);
            free(particles.brightness);
            exit(1);
        }
    

    FILE *file = fopen(filename, "rb");

    if (file == NULL){ // Checks that file got read
        printf("Error reading file");
        exit(1);
    }

    for (int i = 0; i < N; ++i) {
        fread(&particles.x[i], sizeof(double), 1, file);
            fread(&particles.y[i], sizeof(double), 1, file);
            fread(&particles.mass[i], sizeof(double), 1, file);
            fread(&particles.vx[i], sizeof(double), 1, file);
            fread(&particles.vy[i], sizeof(double), 1, file);
            fread(&particles.brightness[i], sizeof(double), 1, file);
    }

    fclose(file);

    print_particles(particles, N);

}