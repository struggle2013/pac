#include <cmath>
#include <cstdio>
#include <omp.h>

struct ParticleSet { 
	float *x, *y, *z;
	float *vx, *vy, *vz; 
};

void MoveParticles(const int nParticles, ParticleSet& particle, const float dt) {
	const int tileSize=16;
	//assert(nParticles%tileSize==0);
	// Loop over particles that experience force
#pragma omp parallel for schedule(guided)
	for (int ii = 0; ii < nParticles; ii+=tileSize) { 

		// Components of the gravity force on particle i
		//float Fx = 0, Fy = 0, Fz = 0; 
		float Fx[tileSize],Fy[tileSize],Fz[tileSize];
		Fx[:]=Fy[:]=Fz[:]=0.0f;
		// Loop over particles that exert force: vectorization expected here

#pragma unroll(tileSize)
		for (int j = 0; j < nParticles; j++) { 

			// Avoid singularity and interaction with self
			const float softening = 1e-20;

#pragma vector aligned
			for(int i=ii;i<ii+tileSize;i++){
				// Newton's law of universal gravity
				const float dx = particle.x[j] - particle.x[i];
				const float dy = particle.y[j] - particle.y[i];
				const float dz = particle.z[j] - particle.z[i];
				const float drSquared  = dx*dx + dy*dy + dz*dz + softening;
				// const float drPower32  = pow(drSquared, 3.0/2.0);
				const float rr = 1.0f/sqrtf(dx*dx+dy*dy+dz*dz+softening);
				const float drPower32 = rr*rr*rr;
				// Calculate the net force
				Fx[i-ii] += dx * drPower32;  
				Fy[i-ii] += dy * drPower32;  
				Fz[i-ii] += dz * drPower32;
			}
		}

		// Accelerate particles in response to the gravitational force
		particle.vx[ii:tileSize] += dt*Fx[0:tileSize]; 
		particle.vy[ii:tileSize] += dt*Fy[0:tileSize]; 
		particle.vz[ii:tileSize] += dt*Fz[0:tileSize];
	}

	// Move particles according to their velocities
#pragma simd
#pragma vector aligned
	for (int i = 0 ; i < nParticles; i++) { 
		particle.x[i]  += particle.vx[i]*dt;
		particle.y[i]  += particle.vy[i]*dt;
		particle.z[i]  += particle.vz[i]*dt;
	}
}

int main(const int argc, const char** argv) {

	// Problem size and other parameters
	const int nParticles = (argc > 1 ? atoi(argv[1]) : 16384);
	const int nSteps = 10;  // Duration of test
	const float dt = 0.01f; // Particle propagation time step

	// Particle data stored as an Array of Structures (AoS)
	// this is good object-oriented programming style,
	// but inefficient for the purposes of vectorization
	//ParticleType* particle = new ParticleType[nParticles];
	ParticleSet particle;
	particle.x=(float*)_mm_malloc(sizeof(float)*nParticles,64);
	particle.y=(float*)_mm_malloc(sizeof(float)*nParticles,64);
	particle.z=(float*)_mm_malloc(sizeof(float)*nParticles,64);
	particle.vx=(float*)_mm_malloc(sizeof(float)*nParticles,64);
	particle.vy=(float*)_mm_malloc(sizeof(float)*nParticles,64);
	particle.vz=(float*)_mm_malloc(sizeof(float)*nParticles,64);
	// Initialize random number generator and particles
	srand(0);
	for(int i = 0; i < nParticles; i++) {
		particle.x[i] = rand()/RAND_MAX;
		particle.y[i] = rand()/RAND_MAX;
		particle.z[i] = rand()/RAND_MAX;
		particle.vx[i] = rand()/RAND_MAX;
		particle.vy[i] = rand()/RAND_MAX;
		particle.vz[i] = rand()/RAND_MAX;
	}

	// Perform benchmark
	int num = omp_get_max_threads();
	printf("\n\033[1mNBODY Version 00\033[0m\n");
	printf("\nPropagating %d particles using %d thread on %s...\n\n", 
			nParticles,num,
#ifndef __MIC__
			"CPU"
#else
			"MIC"
#endif
		  );
	double rate = 0, dRate = 0; // Benchmarking data
	const int skipSteps = 3; // Skip first iteration is warm-up on Xeon Phi coprocessor
	printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);
	for (int step = 1; step <= nSteps; step++) {

		const double tStart = omp_get_wtime(); // Start timing
		MoveParticles(nParticles, particle, dt);
		const double tEnd = omp_get_wtime(); // End timing

		const float HztoInts   = float(nParticles)*float(nParticles-1) ;
		const float HztoGFLOPs = 20.0*1e-9*float(nParticles)*float(nParticles-1);

		if (step > skipSteps) { // Collect statistics
			rate  += HztoGFLOPs/(tEnd - tStart); 
			dRate += HztoGFLOPs*HztoGFLOPs/((tEnd - tStart)*(tEnd-tStart)); 
		}

		printf("%5d %10.3e %10.3e %8.1f %s\n", 
				step, (tEnd-tStart), HztoInts/(tEnd-tStart), HztoGFLOPs/(tEnd-tStart), (step<=skipSteps?"*":""));
		fflush(stdout);
	}
	rate/=(double)(nSteps-skipSteps); 
	dRate=sqrt(dRate/(double)(nSteps-skipSteps)-rate*rate);
	printf("-----------------------------------------------------\n");
	printf("\033[1m%s %4s \033[42m%10.1f +- %.1f GFLOP/s\033[0m\n",
			"Average performance:", "", rate, dRate);
	printf("-----------------------------------------------------\n");
	printf("* - warm-up, not included in average\n\n");
	//delete particle;
	_mm_free(particle.x);
	_mm_free(particle.y);
	_mm_free(particle.z);
	_mm_free(particle.vx);
	_mm_free(particle.vy);
	_mm_free(particle.vz);
}


