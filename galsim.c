#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Introducing the 2D vectorspace
typedef struct{
    double x, y;
} vector2D;

// Distance vector [m]
vector2D r_ij(vector2D u, vector2D v) {
    return (vector2D) {(u.x - v.x), (u.y - v.y)};
}

// Distance scalar [m]
double r_scalar(vector2D u, vector2D v) {
    return sqrt(pow((u.x - v.x),2) + pow((u.y - v.y),2));
}

// Hatfuntion of r [m]
vector2D r_hat(vector2D u, vector2D v) {
    int denum = r_scalar(u,v);
    return (vector2D) {(u.x - v.x)/denum, (u.y - v.y)/denum};
}

// Defining the force finction [N]
vector2D F_ij(int m_i, int m_j, vector2D u, vector2D v, int N) {
    // To mke F_i, make a forloop in main
    double const epsi_0 = pow(10, -3);
    double const G = 100/N;
    return (vector2D){-(G*m_i*m_j * r_ij(u,v).x)/pow(r_scalar(u,v)+epsi_0, 3), -(G*m_i*m_j * r_ij(u,v).y)/pow(r_scalar(u,v)+epsi_0,3)};
}

// Function to print all particles (From ChatGBT)
void print_particles(double particle_info[][6], int N) {
    printf("\nParticle Data (x, y, mass, vx, vy, brightness):\n");
    for (int i = 0; i < N; i++) {
        printf("Particle %d: %.3f %.3f %.3f %.3f %.3f %.3f\n",
               i, particle_info[i][0], particle_info[i][1], particle_info[i][2],
               particle_info[i][3], particle_info[i][4], particle_info[i][5]);
    }
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Wrong input variables");
        return 1;
    }
    // Using converitng the imput to variables
    // argv[0] = ./...
    int const N = atoi(argv[1]);          //The number of stars/particles
    const char *filename = argv[2];       //The filename
    int const nsteps = atoi(argv[3]);     //The number of timesteps
    double const delta_t = atof(argv[4]); //The time step
    int const graphics = atoi(argv[5]);   //On or of, 1 or 0

    // Reeding the file
    double particle_info[N][6];
    FILE *file = fopen(filename, "rb");
    for (int i = 0; i<N; i++) {
        fread(&particle_info[i][0], sizeof(double), 1, file); // Position x
        fread(&particle_info[i][1], sizeof(double), 1, file); // Position y
        fread(&particle_info[i][2], sizeof(double), 1, file); // Mass
        fread(&particle_info[i][3], sizeof(double), 1, file); // Velocity x
        fread(&particle_info[i][4], sizeof(double), 1, file); // Velocity y
        fread(&particle_info[i][5], sizeof(double), 1, file); // Brightness
    }

    print_particles(particle_info, N);

    double F_i[N][2];
    double a_i[N][2];
    vector2D u;
    vector2D v;
    double m_i;
    double m_j;
    for (int i=0; i<N; i++) {
        vector2D F = (vector2D) {0,0}; 
        u = (vector2D) {particle_info[i][0], particle_info[i][1]};
        m_i = particle_info[i][2];
        for (int j=0; j<N; j++) {
            if (i != j) {
                v = (vector2D) {particle_info[j][0], particle_info[j][1]};
                m_j = particle_info[j][2];
                F.x += F_ij(m_i, m_j, u, v, N).x;
                F.y += F_ij(m_i, m_j, u, v, N).y;
            }
        }
        F_i[i][0] = F.x;
        F_i[i][1] = F.y;
        a_i[i][0] = F.x / m_i;
        a_i[i][1] = F.y / m_i;

        // Updating the velocity
        particle_info[i][3] = particle_info[i][3] + delta_t*a_i[i][0];
        particle_info[i][4] = particle_info[i][4] + delta_t*a_i[i][1];

        // Updating the position
        particle_info[i][0] = particle_info[i][0] + delta_t*particle_info[i][3];
        particle_info[i][1] = particle_info[i][1] + delta_t*particle_info[i][4];
    }

    print_particles(particle_info, N);
    
    // Terminal
    // gcc -o galsim galsim.c
    //./galsim 4 input_data/circles_N_4.gal 10 0.001 0

    return 0;
    

}