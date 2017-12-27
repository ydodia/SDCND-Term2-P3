/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

//#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[])
{
	num_particles = 5;
	sig_x = std[0];
	sig_y = std[1];
	sig_theta = std[2];
	default_random_engine gen;
	dist_x = normal_distribution<double>(0, sig_x);
	dist_y = normal_distribution<double>(0, sig_y);
	dist_theta = normal_distribution<double>(0, sig_theta);

	for(int k=0; k<num_particles; ++k)
	{
		Particle p = {};
		p.id = k;
		// generate using gaussian distribution
		p.x = x + dist_x(gen);
		p.y = y + dist_y(gen);
		p.theta = theta + dist_theta(gen);
		p.weight = 1.0;
		particles.push_back(p);
		p.associations.clear();
		p.sense_x.clear();
		p.sense_y.clear();
	}
	weights = vector<double>(num_particles, 1.0);
	is_initialized = true;
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
	if(abs(yaw_rate) > 0.0)
		for(auto& p : particles)
		{
			p.x += velocity / yaw_rate * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta)) + dist_x(gen);
			p.y += velocity / yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate*delta_t)) + dist_y(gen);
			p.theta += yaw_rate * delta_t + dist_theta(gen);
		}
	else
		for(auto& p : particles)
		{
			p.x += velocity * delta_t * cos(p.theta) + dist_x(gen);
			p.y += velocity * delta_t * sin(p.theta) + dist_y(gen);
			p.theta += dist_theta(gen);
		}
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
}

void ParticleFilter::dataAssociation(const std::vector<LandmarkObs>& predicted, const Map& map_landmarks, double sensor_range)
{

	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map& map_landmarks)
{
	static double denom = 2.0 * M_PI * sig_x * sig_y;
	for(auto& p : particles)
	{
		// 1 - Find predicted set of observations that are within sensor range
		vector<LandmarkObs> predicted;
		for(const auto& obs : observations)
		{
			// convert from vehicle to map coords
			double x_map = (obs.x * cos(p.theta)) - (sin(p.theta) * obs.y);
			double y_map = (obs.x * sin(p.theta)) + (cos(p.theta) * obs.y);
			double range = sqrt( pow(x_map, 2.0) + pow(y_map, 2.0) );
			if (range <= sensor_range)
				predicted.push_back(LandmarkObs{obs.id, x_map + p.x, y_map + p.y});
		}
		// 2 - If we have equal num of predicted to measured, we then associate to closest landmark
		if(predicted.size() == observations.size())
		{
			for(const auto& obs : predicted)
			{
				double min_dist = sensor_range;
				int idx = 0;
				double x_pos{0}, y_pos{0};
				for(const auto& landmark : map_landmarks.landmark_list)
				{
					double x = obs.x - landmark.x_f;
					double y = obs.y - landmark.y_f;
					double dist = sqrt(x*x + y*y);
					if(dist < min_dist)
					{
						min_dist = dist;
						idx = landmark.id_i;
						x_pos = landmark.x_f;
						y_pos = landmark.y_f;
					}
				}
				// 3 - Calculate weight for each association and then total.
				double pwr = pow(obs.x - x_pos, 2.0)  / (2.0 * pow(sig_x, 2.0));
				pwr += pow(obs.y - y_pos, 2.0)  / (2.0 * pow(sig_y, 2.0));
				p.weight *= exp(-1.0 * pwr) / denom;
			}
		}
		else
			p.weight = 0;
	}

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
}

void ParticleFilter::resample()
{
	for(int k=0; k<num_particles; ++k)
		weights[k] = particles[k].weight;
	dist_discrete = discrete_distribution<int>(weights.begin(), weights.end());

	vector<Particle> resampled_particles;
	for(int k=0; k<num_particles; ++k)
		resampled_particles.push_back( particles[dist_discrete(gen)] );
	particles.swap(resampled_particles);

	for(auto& p : particles)
		p.weight = 1.0;


	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
