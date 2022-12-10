//This is mainly based on the sequence given by https://conference.sdo.esoc.esa.int/proceedings/sdc7/paper/14/SDC7-paper14.pdf
//The timestep calculation is based on Gravitational N-Body Simulations, SVERRE J. AARSETH

#include <iostream>
#include <map>
#include <vector>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <exception>
#include <sstream>
#include <time.h>
#include "hrtime.h"
#include "omp.h"
using namespace std;

const long double G = 6.6743015e-11, PI = 3.1415926535898;
vector <long double> add_init_position, add_init_velocity, burn_vector, burn_ori, burn_ori_rate;
vector <string> load_results;
long double add_mass;
//long double r_dot_v_relative, abs_r_rel, a_comp, a_dot_comp_1, a_dot_comp_2, abs_v_rel, comp_a, comp_b, comp_c, body_timestep, abs_a_t_dot, abs_a_d_dot, abs_a_dot, abs_a;
string line, add_id, acting_id, burn_id, burn_body;
//int i, itts, 
int time_counter;
int burn_count = 0;
//pair<string, vector<long double> > burn_values, burn_values_next;


class System {
public:
	map<string, vector<long double> > bodies;
	//map<string, vector<long double> > bodies_next;
	map<string, pair<string, vector<long double> > > burns;
	vector<string> body_ids;
	vector<vector<long double>> body_info, body_info_next;
	long double accuracy, timestep, time, next_timestep;
	int output_rate;
	string info, debug_info;

	void SetAccuracy(int new_accuracy) {
		accuracy = new_accuracy;
	}

	void SetStartTime(int start_time) {
		time = start_time;
	}

	void SetOutToFile(const char* file_name, int rate) {
		#pragma warning(suppress : 4996)
		freopen(file_name, "w", stdout);
		output_rate = rate;
	}

	void AddBody(string id, long double mass, vector<long double> init_position, vector<long double> init_velocity) {
		bodies[id].clear();//so adding bodies with the same id twice doesn't break it
		bodies[id].push_back(mass);
		int i;
		for (i = 0;i < init_position.size();i++) {
			bodies[id].push_back(init_position.at(i));
		}

		for (i = 0;i < init_velocity.size();i++) {
			bodies[id].push_back(init_velocity.at(i));
		}
	}

	void AddBurn(string burn_id, string body_id, long double start_time, long double end_time, long double acceleration, vector<long double> orientation, vector<long double> orientation_rate) {
		//Gives burn_id:(body_id:[start,end,accel,orie_x,orie_y,orie_z,orie_ra_x...])
		int i;
		burn_vector.clear();
		burn_vector.push_back(start_time);
		burn_vector.push_back(end_time);
		burn_vector.push_back(acceleration);
		for (i = 0;i < 3;i++) {
			burn_vector.push_back(orientation.at(i));
		}
		for (i = 0;i < 3;i++) {
			burn_vector.push_back(orientation_rate.at(i));
		}
		burns[burn_id] = make_pair(body_id, burn_vector);
	}

	void LoadFile(string filename) {
		//Format is #body,id,mass,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z
		int i;
		ifstream file;
		file.open(filename);
		if (file.is_open()) {
			string line;
			while (getline(file, line)) {
				if (line.compare(0, 1, "#") == 0) {
					if (line.compare(1, 4, "Body") == 0) {
						load_results.clear();
						stringstream s_stream(line);
						while (s_stream.good()) {
							string substr;
							getline(s_stream, substr, ',');
							load_results.push_back(substr);
						}
						add_init_position.clear();
						add_init_velocity.clear();

						add_id = load_results.at(1);
						add_mass = stold(load_results.at(2));

						for (i = 0;i < 3;i++) {
							add_init_position.push_back(stold(load_results.at(3 + i)));
							add_init_velocity.push_back(stold(load_results.at(6 + i)));
						}
						AddBody(add_id, add_mass, add_init_position, add_init_velocity);

					}
					if (line.compare(1, 4, "Burn") == 0) {
						load_results.clear();
						stringstream s_stream(line);
						while (s_stream.good()) {
							string substr;
							getline(s_stream, substr, ',');
							load_results.push_back(substr);
						}
						burn_ori.clear();burn_ori_rate.clear();
						for (i = 0;i < 3;i++) {
							burn_ori.push_back(stold(load_results.at(5 + i)));
							burn_ori_rate.push_back(stold(load_results.at(8 + i)));
						}
						burn_id = load_results.at(1);
						AddBurn(to_string(burn_count), load_results.at(1), stold(load_results.at(2)), stold(load_results.at(3)), stold(load_results.at(4)), burn_ori, burn_ori_rate);
						burn_count++;
					}
					if (line.compare(1, 8, "Addition") == 0) {
						cout << "This functionality has not yet been added";
						throw exception();
					}
				}
			}
			file.close();
		}

		else {
			cout << "cant not open!";
		}
	}

	void initialize() {
		body_ids.clear();
		body_info.clear();
		for (std::pair<std::string, vector <long double> > body_itterator : bodies) {
			body_ids.push_back(body_itterator.first);
			body_info.push_back(body_itterator.second);
		}
		body_info_next = body_info;

	}

	//vector<long double> CalcAAndDot(string body_id, vector<long double> r, vector<long double> v) {
	//	a.clear();
	//	a_dot.clear();
	//	for (i = 0;i < 3;i++) {
	//		a.push_back(0);
	//		a_dot.push_back(0);
	//	}
	//	for (std::pair<std::string, vector <long double> > body_itterator : bodies) {
	//		acting_id = body_itterator.first;
	//		if (acting_id != body_id) {
	//			acting_values = body_itterator.second;

	//			r_rel.clear();
	//			v_rel.clear();
	//			for (i = 0;i < 3;i++) {
	//				r_rel.push_back(r.at(i) - acting_values.at(1 + i));
	//			}
	//			for (i = 0;i < 3;i++) {
	//				v_rel.push_back(v.at(i) - acting_values.at(4 + i));
	//			}
	//			for (i = 0;i < 3;i++) {
	//				r_dot_v_relative = r_rel.at(i) * v_rel.at(i);
	//			}

	//			abs_r_rel = sqrt(pow(r_rel.at(0), 2) + pow(r_rel.at(1), 2) + pow(r_rel.at(2), 2));

	//			a_comp = -G * acting_values.at(0) / pow(abs_r_rel, 3);
	//			a_dot_comp_1 = 3 * G * acting_values.at(0) * r_dot_v_relative / pow(abs_r_rel, 5);
	//			a_dot_comp_2 = -G * acting_values.at(0) / pow(abs_r_rel, 3);

	//			for (i = 0;i < 3;i++) {
	//				a.at(i) = a.at(i) + a_comp * r_rel.at(i);
	//				a_dot.at(i) = a_dot.at(i) + a_dot_comp_1 * r_rel.at(i) + a_dot_comp_2 * v_rel.at(i);
	//			}
	//		}
	//	}
	//	//TODO: Swap round looking at body id first with looking at time
	//	if (burns.size() > 0) {
	//		for (pair<string, pair<string, vector<long double> > > burn_itt : burns) {
	//			burn_id = burn_itt.first;
	//			burn_values = burn_itt.second;
	//			burn_body = burn_values.first;
	//			if (burn_body == body_id) {
	//				burn_vector = burn_values.second;
	//				if (time > burn_vector.at(0) && time < burn_vector.at(1)) {
	//					for (i = 0;i < 3;i++) {
	//						a.at(i) = burn_vector.at(2) * burn_vector.at(3 + i);
	//					}
	//					burn_vector_next.clear();
	//					for (i = 0;i < 3;i++) {
	//						burn_vector_next.push_back(burn_vector.at(i));
	//					}
	//					for (i = 0;i < 3;i++) {
	//						burn_vector_next.push_back(burn_vector.at(3 + i) + burn_vector.at(6 + i) * timestep);
	//					}
	//					for (i = 0;i < 3;i++) {
	//						burn_vector_next.push_back(burn_vector.at(6 + i));
	//					}
	//				}
	//			}
	//		}
	//	}
	//	
	//	
	//	for (i = 0; i < 3; i++) {
	//		oupt
	//	}
	//		
	//}


	void CalAAndADot(long double* output_vec, string id, long double* r, long double* v) {
		long double a[3], a_dot[3], r_rel[3], v_rel[3];
		long double r_dot_v_relative, abs_r_rel, a_comp, a_dot_comp_1, a_dot_comp_2;
		string acting_id;
		vector <long double> acting_values;
		int i;

		for (i = 0; i < 3; i++) {
			a[i] = 0;
			a_dot[i] = 0;
		}

		for (i = 0; i < body_ids.size(); i++) {
			acting_id = body_ids[i];
			if (acting_id != id) {
				acting_values = body_info[i];

				int j;
				for (j = 0; j < 3; j++)
					r_rel[j] = r[j] - acting_values[1+j];
				for (j = 0; j < 3; j++)
					v_rel[j] = v[j] - acting_values[4+j];
				for (j = 0; j < 3; j++)
					r_dot_v_relative = r_rel[j] * v_rel[j];

				abs_r_rel = sqrt(pow(r_rel[0], 2) + pow(r_rel[1], 2) + pow(r_rel[2], 2));

				a_comp = -G * acting_values[0] / pow(abs_r_rel, 3);
				a_dot_comp_1 = 3 * G * acting_values[0] * r_dot_v_relative / pow(abs_r_rel, 5);
				a_dot_comp_2 = -G * acting_values[0] / pow(abs_r_rel, 3);

				for (j = 0; j < 3; j++) {
					a[j] = a[j] + a_comp * r_rel[j];
					a_dot[j] = a_dot[j] + a_dot_comp_1 * r_rel[j] + a_dot_comp_2 * v_rel[j];
				}

			}
		}

		if (burns.size() > 0) {
			for (pair<string, pair<string, vector<long double> > > burn_itt : burns) {
				string burn_body = burn_itt.second.first;
				
				if (burn_body == id) {
					vector <long double> burn_value = burn_itt.second.second;
					if (time > burn_value.at(0) && time < burn_value.at(1)) {
						for (i = 0;i < 3;i++) {
							a[i] = burn_value.at(2) * burn_value.at(3 + i);
						}
						/*burn_vector_next.clear();
						for (i = 0;i < 3;i++) {
							burn_vector_next.push_back(burn_vector.at(i));
						}
						for (i = 0;i < 3;i++) {
							burn_vector_next.push_back(burn_vector.at(3 + i) + burn_vector.at(6 + i) * timestep);
						}
						for (i = 0;i < 3;i++) {
							burn_vector_next.push_back(burn_vector.at(6 + i));
						}*/
					}
				}
			}
		}

		for (i = 0; i < 3; i++) {
			output_vec[i] = a[i];
			output_vec[3 + i] = a_dot[i];
		}

	}
	

	void Step() {
		next_timestep = 9999999999;
		#pragma omp parallel for 
		for (int i = 0; i < body_ids.size(); i++) {

			int j, itts;
			long double r_0[3], v_0[3];
			for (j = 0; j < 3; j++) {
				r_0[j] = body_info[i][1+j];
				v_0[j] = body_info[i][4+j];
			}

			string body_id = body_ids[i];
			long double output_vec[6];
			CalAAndADot(output_vec, body_id, r_0, v_0);

			long double a_0[3], a_0_dot[3];
			for (j = 0; j < 3; j++) {
				a_0[j] = output_vec[j];
				a_0_dot[j] = output_vec[j + 3];

			}

			long double r_p[3], v_p[3];
			for (j = 0; j < 3; j++) {
				r_p[j] = r_0[j] + v_0[j] * timestep + 0.5 * a_0[j] * pow(timestep, 2) + (1 / 6) * a_0_dot[j] * pow(timestep, 3);
				v_p[j] = v_0[j] + a_0[j] * timestep + 0.5 * a_0_dot[j] * pow(timestep, 2);
			}


			long double a_p[3], a_p_dot[3];
			for (itts = 0;itts < 2;itts++) {
				CalAAndADot(output_vec, body_id, r_p, v_p);
				for (j = 0; j < 3;j++) {
					a_p[j] = output_vec[j];
					a_p_dot[j] = output_vec[j + 3];
				}

				for (j = 0; j < 3; j++) {
					v_p[j] = v_0[j] + 0.5 * (a_0[j] + a_p[j]) * timestep + (1 / 12) * (a_0_dot[j] - a_p_dot[j] * pow(timestep, 2));
					r_p[j] = r_0[j] + 0.5 * (v_p[j] + v_0[j]) * timestep + (1 / 12) * (a_0[j] - a_p[j]) * pow(timestep, 2);
				}
			}


			vector<long double> body_final;
			body_final.push_back(body_info[i][0]);
			for (j = 0; j < 3; j++) {
				body_final.push_back(r_p[j]);
			}
			for (j = 0; j < 3; j++) {
				body_final.push_back(v_p[j]);
			}

			body_info_next[i] = body_final;


		}
		body_info = body_info_next;
		time += timestep;
	}

	void Output() {
		if (time_counter == output_rate) {
			cout << "#" + to_string(time) + "\n";
			string body_id;
			vector<long double> body_values;
			for (int i = 0; i < body_ids.size(); i++) {
				body_id = body_ids[i];
				body_values = body_info[i];
				cout<<body_id;
				for (int j = 0; j < body_values.size(); j++) {
					cout << "," + to_string(body_values.at(j));
				}
				cout << "\n";
			}
			time_counter = 0;
		}
		else {
			time_counter++;
		}
	}
};



int main(int argc, char * argv[]) {
	int num_threads = 8;

	if(argc > 2) {
		std::cout << "Usage:" << argv[0] << " [num_threads]\n";
	}

	if(argc == 2) num_threads = atoi(argv[1]);

	omp_set_dynamic(0);
	omp_set_num_threads(num_threads);


	System sys;
	sys.SetStartTime(0);
	sys.timestep = 0.25;
	sys.LoadFile("./test/solar_system.start");
	sys.SetOutToFile("./test/test_OpenMP.txt", 10);
	sys.initialize();

	double endtime = 30000;
	double steps = endtime / sys.timestep;
 

	double elapsed = 0;
	
	double start, end;
	for (int step = 0; step < steps; step++) {
	sys.Output();
	start = getElapsedTime();
	sys.Step();
	end = getElapsedTime(); 
	elapsed += (end - start);
	}

	fprintf(stdout, "\n==============================\nCode took: %f seconds\n", elapsed);

	fclose(stdout);
	stdout = fdopen(0, "w");
	fprintf(stdout, "Code took: %f seconds\n", elapsed);

}
