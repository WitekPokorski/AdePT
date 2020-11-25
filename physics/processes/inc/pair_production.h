#ifndef PAIRPRODUCTIONH
#define PAIRPRODUCTIONH

#include "process.h"

class pair_production: public process  {
    public:
        __device__ pair_production() {}
        __device__ virtual float GetPhysicsInteractionLength(int particle_index, adept::BlockData<track> *block, curandState_t *states) const;
        __device__ virtual void GenerateInteraction(int particle_index, adept::BlockData<track> *block, curandState_t *states);
};

__device__ float pair_production::GetPhysicsInteractionLength(int particle_index, adept::BlockData<track> *block, curandState_t *states) const {
    
    float current_length = curand_uniform(states) * 100.0f; // here I need to calculate the IL based on the particle energy, material, etc
    return current_length;
}

__device__ void pair_production::GenerateInteraction(int particle_index, adept::BlockData<track> *block, curandState_t * states)
{
    track* mytrack = &((*block)[particle_index]);

    if (mytrack->energy > 0.001f)
    {
        // pair production
        float eloss = 0.5f * mytrack->energy;

        mytrack->energy -= eloss;

        auto secondary_track = block->NextElement();
        assert(secondary_track != nullptr && "No slot available for secondary track");
        secondary_track->energy = eloss;
    }
    else
    {
        mytrack->status = dead;
    }
    

}



#endif