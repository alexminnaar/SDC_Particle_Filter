/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <float.h>
#include <climits>
#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

// double ParticleFilter::distance(double x_1, double y_1, double x_2, double y_2){

//     return sqrt(pow(x_2-x_1,2) + pow(y_2-y_1,2));

// }

void ParticleFilter::init(double x, double y, double theta, double std[])
{
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    default_random_engine gen;

    double std_x = std[0];
    double std_y = std[1];
    double std_psi = std[2];

    //set number of particles
    num_particles = 101;

    //create noise distributions - zero mean gaussian
    normal_distribution<double> dist_x(0, std_x);
    normal_distribution<double> dist_y(0, std_x);
    normal_distribution<double> dist_theta(0, std_x);

    //create particles
    for (int i = 0; i < num_particles; i++)
    {
        Particle rand_particle;

        //set x, y, and theta with some added gaussian noise
        rand_particle.id = i;
        rand_particle.x = x+dist_x(gen);
        rand_particle.y = y+dist_y(gen);
        rand_particle.theta = theta+dist_theta(gen);
        rand_particle.weight = 1.0;

        //add to particle vector
        particles.push_back(rand_particle);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;

    //create noise distributions
    normal_distribution<double> noise_x(0,std_pos[0]);
    normal_distribution<double> noise_y(0,std_pos[1]);
    normal_distribution<double> noise_theta(0,std_pos[2]);

    //I assume we are supposed to fill in the sense fields of the particle
    for (int i = 0; i < num_particles; ++i)
    {
        double x = particles[i].x;
        double y = particles[i].y;
        double theta = particles[i].theta;

        //predict x and y positions based on motion model
        //TODO add gaussian noise

        if(fabs(yaw_rate) < 0.00001){
            particles[i].x += velocity*delta_t*cos(theta);
            particles[i].y += velocity*delta_t*sin(theta);
        }
        else{
            particles[i].x += (velocity / yaw_rate) * (sin(theta + yaw_rate * delta_t) - sin(theta));
            particles[i].y+=(velocity / yaw_rate) * (cos(theta) - cos(theta + yaw_rate * delta_t));
            particles[i].theta += yaw_rate*delta_t;;
        }
        
        //Add some noise
        particles[i].x +=  noise_x(gen);
        particles[i].y +=  noise_y(gen);
        particles[i].theta+= noise_theta(gen);
    }
}



void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations)
{
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    //for each observation, find the closest landmark
    for(vector<LandmarkObs>::iterator obs_it = observations.begin(); obs_it != observations.end(); ++obs_it){
        
        LandmarkObs obs = *obs_it;

        //keep track of closest distance, start with large number
        double closest_dist = numeric_limits<double>::max();
        int closest_id = -1;

        for(vector<LandmarkObs>::iterator pred_it = predicted.begin(); pred_it != predicted.end(); ++pred_it){

            LandmarkObs pred = *pred_it;
            double distance = dist(obs.x,obs.y,pred.x,pred.y);

            //keep track of closest landmark
            if(distance<closest_dist){
                closest_dist=distance;
                closest_id=pred.id;
            }
        }

        //assign observation id to closest landmark
        obs_it->id = closest_id;
        cout<<obs_it->id<<endl;

    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks)
{
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html

    //iterate through particles
    for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it)
    {
        Particle p = *it;

        vector<LandmarkObs> landmarks_in_range;

        //get landmarks within the range of this sensor
        for(int i =0; i< map_landmarks.landmark_list.size();i++){

            int landmark_id = map_landmarks.landmark_list[i].id_i;
            float landmark_x = map_landmarks.landmark_list[i].x_f;
            float landmark_y = map_landmarks.landmark_list[i].x_f;

            //check if the landmark is within the range for this particle
            if(fabs(p.x - landmark_x) <= sensor_range && fabs(p.y - landmark_y)<=sensor_range){
                landmarks_in_range.push_back(LandmarkObs{landmark_id,landmark_x,landmark_y});
            }
        }

        //transform vehicle observations into map coordinates

        //record transformed observations in a vector
        vector<LandmarkObs> transformed_obs;

        //iterate through vehicle observations and transform them into map coordinates
        for(vector<LandmarkObs>::iterator obs_it = observations.begin(); obs_it!=observations.end();++obs_it ){

            //get x and y observations in vehicle coordinates
            double obs_x_vehicle = obs_it->x;
            double obs_y_vehicle = obs_it->y;
            
            //convert observations to map coordinates wrt this particle
            double x_map = obs_x_vehicle*cos(p.theta) - obs_y_vehicle*sin(p.theta)+ p.x;
            double y_map = obs_x_vehicle*sin(p.theta) + obs_y_vehicle*cos(p.theta)+p.y;
            
            //add to transformed vector
            transformed_obs.push_back(LandmarkObs{obs_it->id, x_map, y_map});
        }

        //use data association to find closest landmark for each observation
        dataAssociation(landmarks_in_range,transformed_obs);


        //for each observation/landmark pair, calculate weight for this particle 
        //initially set weight to 1 because we are about to iteratively multiply to get the final weight
        it->weight = 1.0;

        for(vector<LandmarkObs>::iterator trans_obs_it = transformed_obs.begin(); trans_obs_it!= transformed_obs.end();++trans_obs_it){

            LandmarkObs obs = *trans_obs_it;

            double landmark_associated_x, landmark_associated_y;

            //get the associated landmark for this observation which was computed in dataAssociation()
            for(vector<LandmarkObs>::iterator landmark_it = landmarks_in_range.begin(); landmark_it!=landmarks_in_range.end();++landmark_it){

                LandmarkObs landmark = *landmark_it;

                //check if this observation was associated with this landmark
                if(obs.id == landmark.id){
                    landmark_associated_x=landmark.x;
                    landmark_associated_y=landmark.y;
                }
            }

            //add this landmarks contribution to the particle's weight
            double weight_contribution = (1/(2*M_PI*std_landmark[0]*std_landmark[1]))*exp(-(pow(landmark_associated_x - obs.x,2)/(2*pow(std_landmark[0],2)) + (pow(landmark_associated_y - obs.y,2)/(2*pow(std_landmark[1],2))))  );
            //cout<<it->weight;
            it->weight *= weight_contribution;

        }
    
    }
}

void ParticleFilter::resample()
{
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    default_random_engine gen;

    //first create a vector of particles for the next iteration
    vector<Particle> next_particles;

    vector<double> current_weights;

    for(vector<Particle>::iterator p_it= particles.begin(); p_it!=particles.end();++p_it){
        current_weights.push_back(p_it->weight);
    }

    //create a distribution over particle indexes
    uniform_int_distribution<int> particle_index_dist(0, num_particles-1);
    int particle_index = particle_index_dist(gen);

    //create a distribution over weight values
    double max_weight = *max_element(current_weights.begin(), current_weights.end());
    uniform_real_distribution<double> weight_dist(0.0, max_weight);

    double beta=0.0;

    //perform weight resampling
    for(int i=0; i < num_particles; i++){
        beta+=2.0*weight_dist(gen);

        while(beta>current_weights[particle_index]){
            beta -= current_weights[particle_index];
            particle_index= (particle_index+1)%num_particles;
        }
        next_particles.push_back(particles[particle_index]);
    }

    //overwrite old particles
    particles = next_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1); // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
    vector<double> v = best.sense_x;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1); // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
    vector<double> v = best.sense_y;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1); // get rid of the trailing space
    return s;
}
