#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

double laneToD(int lane)
{
	const double lane_width_m = 4;
	double d = (lane + 0.5) * lane_width_m; //calculate d, 0.5 is to move to middle of the lane
	return d;
}

double dToLane(double d)
{
	for (int lane = 0; lane < 3; lane++)
	{
		if ((d > laneToD(lane) - 2) && (d < laneToD(lane) + 2))
		{
			return lane;
		}

	}

	return 0; //if we got here we got problems
}

enum LaneDirection
{
	LaneDirection_left,
	LaneDirection_right
};

enum Decisions
{
	Decisions_straight,
	Decisions_left,
	Decisions_right,
};

struct sensorFustionSummary
{
	double leadCarDist = 999;
	double lefttLeadCarDist = 999;
	double leftTailCarDist = 999;
	double rightLeadCarDist = 999;
	double rightTailCarDist = 999;
	double laneSpeed[3] = { 0, 0, 0 };
} SensorFusionSummary;

struct decisionResult
{
	double speed;
	int lane;
} DecisionResult;


int getNextLane(int lane, LaneDirection laneDirection)
{
	int nextLane;

	if (lane == 0)
	{
		if (laneDirection == LaneDirection_left)
		{
			nextLane = -1;
		}
		else
		{
			nextLane = 1;
		}

	}

	if (lane == 1)
	{
		if (laneDirection == LaneDirection_left)
		{
			nextLane = 0;
		}
		else
		{
			nextLane = 2;
		}

	}

	if (lane == 2)
	{
		if (laneDirection == LaneDirection_left)
		{
			nextLane = 1;
		}
		else
		{
			nextLane = -1;
		}

	}

	return nextLane;
}

sensorFustionSummary summarizeSensorFusion(int lane, double car_s, json sensor_fusion, double prev_size, double dt_s)
{
	sensorFustionSummary result = sensorFustionSummary();

	//search through sensor fusion vector
	for (int i = 0; i < sensor_fusion.size(); i++)
	{

		//pulls some info from the car in sensor fusion
		double vx = sensor_fusion[i][3];
		double vy = sensor_fusion[i][4];
		double check_speed = sqrt(pow(vx, 2) + pow(vy, 2));
		double check_car_s = sensor_fusion[i][5];
		double check_car_s_future = check_car_s + prev_size* dt_s * check_speed;	//extrapolate where the car will be in the furture
		float d = sensor_fusion[i][6];
		double check_lane = dToLane((double)d);
		double check_distance = check_car_s_future - car_s;


		//Define left and right lanes, if invalid, they are set to -1
		int leftLane = getNextLane(lane, LaneDirection_left);
		int rightLane = getNextLane(lane, LaneDirection_right);

		if (check_lane == lane)
		{
			//if in the lane check distance to lead car
			if (check_distance > 0)
			{
				//find the closest car and log distance and speed

				if (abs(result.leadCarDist) > check_distance)
				{
					result.leadCarDist = abs(check_distance);
					result.laneSpeed[lane] = check_speed;
				}
				
			}
			
		}
		else if (check_lane == leftLane)
		{
			//if in the lane check distance to lead car
			if (check_distance > 0)
			{
				//find the closest car and log distance and speed

				if (abs(result.leadCarDist) > check_distance)
				{
					result.leftTailCarDist = abs(check_distance);
					result.laneSpeed[leftLane] = check_speed;
				}

			}

			//if in the lane check distance to tail car
			if (check_distance <= 0)
			{
				//find the closest car and log distance and speed

				if (abs(result.leadCarDist) > check_distance)
				{
					result.leftTailCarDist = abs(check_distance);
					result.laneSpeed[leftLane] = check_speed;
				}

			}
		}
		else if (check_lane == rightLane)
		{
			//if in the lane check distance to lead car
			if (check_distance > 0)
			{
				//find the closest car and log distance and speed

				if (abs(result.leadCarDist) > check_distance)
				{
					result.rightLeadCarDist = abs(check_distance);
					result.laneSpeed[rightLane] = check_speed;
				}

			}

			//if in the lane check distance to tail car
			if (check_distance <= 0)
			{
				//find the closest car and log distance and speed

				if (abs(result.leadCarDist) > check_distance)
				{
					result.rightTailCarDist = abs(check_distance);
					result.laneSpeed[rightLane] = check_speed;
				}

			}
		}
	}

	return result;
}

decisionResult Decide(sensorFustionSummary summary, int lane, double ref_vel_mph, double dt_s)
{
	//define thresholds
	double sameLaneDist_m = 30;
	double laneChaneDist_m = 30;
	double openLaneDist_m = 60;
	double target_v = 49.5;
	double straightBonus = 0;


	//Check for collisions when merging left or right

	bool allowMergeLeft = true;
	bool allowMergeRight = true;

	if ((summary.lefttLeadCarDist < laneChaneDist_m)||(getNextLane(lane, LaneDirection_left)<0))
	{
		allowMergeLeft = false;
	}
	if ((summary.leftTailCarDist < laneChaneDist_m)||(getNextLane(lane, LaneDirection_left)<0))
	{
		allowMergeLeft = false;
	}

	if ((summary.rightLeadCarDist< laneChaneDist_m)||(getNextLane(lane, LaneDirection_right)<0))
	{
		allowMergeRight = false;
	}
	if ((summary.rightTailCarDist < laneChaneDist_m)||(getNextLane(lane, LaneDirection_right)<0))
	{
		allowMergeRight = false;
	}

	//calculate speed cost
	double speedCost[] = { 0,0,0 };

	for (int i = 0; i < 3; i++)
	{
		speedCost[i] = abs(summary.laneSpeed[i] - target_v);
	}


	//calculate open lane bonus
	double openBonusStraight = summary.leadCarDist - openLaneDist_m;
	if (openBonusStraight < 0)
	{
		openBonusStraight = 0;
	}

	double openBonusLeft = summary.lefttLeadCarDist - openLaneDist_m;
	if (openBonusLeft < 0)
	{
		openBonusLeft = 0;
	}

	double openBonusRight = summary.rightLeadCarDist - openLaneDist_m;
	if (openBonusRight < 0)
	{
		openBonusRight = 0;
	}


	double straightCost = speedCost[lane] - straightBonus- openBonusStraight;
	double leftCost = speedCost[getNextLane(lane, LaneDirection_left)]- openBonusLeft;
	double rightCost = speedCost[getNextLane(lane, LaneDirection_right)] - openBonusRight;

	


	//start with going straight
	Decisions bestDecision = Decisions_straight;
	double bestCost = straightCost;
	int bestLane = lane;

	//check left
	if (allowMergeLeft)
	{
		if (leftCost < bestCost)
		{
			bestDecision = Decisions_left;
			bestCost = leftCost;
			bestLane = getNextLane(lane, LaneDirection_left);
		}
	}

	//check right
	if (allowMergeRight)
	{
		if (rightCost < bestCost)
		{
			bestDecision = Decisions_right;
			bestCost = rightCost;
			bestLane = getNextLane(lane, LaneDirection_right);
		}
	}


	//if we haven't gotten too close, then speed back up to reference velocity
	double max_accel_mpss = 9.5;
	double mphTomps = 0.447;
	double new_speed = ref_vel_mph;
	int new_lane = lane;



	if (summary.leadCarDist < sameLaneDist_m)
	{
		new_lane = bestLane;

		//slow down if we got too close to the car in our lane
		if (bestDecision == Decisions_straight)
		{
			//slow down if our best decision is to stay in the same lane
			new_speed = ref_vel_mph - max_accel_mpss * dt_s / mphTomps;
			if (ref_vel_mph < summary.laneSpeed[lane]) // follow at slightly less than the leading car's speed
			{
				new_speed = summary.laneSpeed[lane] * 0.9;
			}
		}
	}
	else if (ref_vel_mph < target_v)
	{
		new_speed = ref_vel_mph + max_accel_mpss * dt_s / mphTomps;
	}


	decisionResult result = decisionResult();
	result.lane = new_lane;
	result.speed = new_speed;

	


	return result;
	

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  int lane = 1;
  double ref_vel_mph = 0;

#ifdef UWS_VCPKG
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &lane, &ref_vel_mph](uWS::WebSocket<uWS::SERVER> *ws, char *data, size_t length,
                     uWS::OpCode opCode) {
#else
  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy, &lane, &ref_vel_mph](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
	  uWS::OpCode opCode) {
#endif
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;

	double dt_s = 0.02;
	double mphTomps = 0.447;

    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];


          	json msgJson;

          	

			double ref_x = car_x;
			double ref_y = car_y;
			double ref_yaw = deg2rad(car_yaw);
			int prev_size = previous_path_x.size();


			//////////////////////////////////////////////////////////////////////////////////////////////////
			//sensor fusion

			//take the s value from the end of the path
			if (prev_size > 0)
			{
				car_s = end_path_s;
			}

			sensorFustionSummary summary =  summarizeSensorFusion(lane, car_s, sensor_fusion, prev_size, dt_s);
			decisionResult decision = Decide(summary, lane, ref_vel_mph, dt_s);

			lane = decision.lane;
			ref_vel_mph = decision.speed;

			//cout << lane << endl;
			//cout << ref_vel_mph << endl;
			//cout << "-----------" << endl;


			///////////////////////////////////////////////////////////////////////////////////////////////////
			//Build spline waypoints

			vector<double> spline_ptsx;
			vector<double> spline_ptsy;

			
			if (prev_size < 2) //if almost empty
			{
				//estimate previous location
				double prev_car_x = car_x - cos(car_yaw);
				double prev_car_y = car_y - sin(car_yaw);

				//add the guess of where we have been
				spline_ptsx.push_back(prev_car_x);
				spline_ptsx.push_back(car_x);

				spline_ptsy.push_back(prev_car_y);
				spline_ptsy.push_back(car_y);
			}
			else
			{
				//calculate previous location
				ref_x = previous_path_x[prev_size - 1];
				ref_y = previous_path_y[prev_size - 1];

				double prev_ref_x = previous_path_x[prev_size - 2];
				double prev_ref_y = previous_path_y[prev_size - 2];
				
				//calculate yaw
				ref_yaw = atan2(ref_y - prev_ref_y, ref_x - prev_ref_x);

				//add the end of the path
				spline_ptsx.push_back(prev_ref_x);
				spline_ptsx.push_back(ref_x);

				spline_ptsy.push_back(prev_ref_y);
				spline_ptsy.push_back(ref_y);
			}

			//Add sparse points along the path
			int spline_numSteps = 3;
			int spline_step_m = 30;

			for (int i = 0; i < spline_numSteps; i++)
			{
				double next_s = car_s + spline_step_m * (i + 1);
				double next_d = laneToD(lane);

				vector<double> next_xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

				double next_x = next_xy[0];
				double next_y = next_xy[1];

				spline_ptsx.push_back(next_x);
				spline_ptsy.push_back(next_y);
			}

			//bring the points into the car's reference frame
			for (int i = 0; i < spline_ptsx.size(); i++)
			{
				double shift_x = spline_ptsx[i] - ref_x;
				double shift_y = spline_ptsy[i] - ref_y;

				spline_ptsx[i] = (shift_x*cos(-ref_yaw) - shift_y * sin(-ref_yaw));
				spline_ptsy[i] = (shift_x*sin(-ref_yaw) + shift_y * cos(-ref_yaw));
			}

			////////////////////////////////////////////////////
			//create the spline
			tk::spline s;
			s.set_points(spline_ptsx, spline_ptsy);

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////

			vector<double> next_x_vals;
			vector<double> next_y_vals;

			int path_size = previous_path_x.size();

			//add the previous path points
			for (int i = 0; i < path_size; i++)
			{
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}

			//calculate the spacing along the spline for target velocity
			
			double horizon_x = 30;
			double horizon_y = s(horizon_x);
			double target_dist = sqrt(pow(horizon_x, 2) + pow(horizon_y, 2)); //pythagorean theorem, estimate the path as a line for this
			double num_points = target_dist / (dt_s*ref_vel_mph*mphTomps);

			double carRef_x = 0;
			int numNextPts = 50;

			for (int i = 0; i < numNextPts - path_size; i++)
			{
				carRef_x += horizon_x / num_points;
				double carRef_y = s(carRef_x);

				//transform back into map space
				double rotate_x = (carRef_x*cos(ref_yaw) - carRef_y * sin(ref_yaw));
				double rotate_y = (carRef_x*sin(ref_yaw) + carRef_y * cos(ref_yaw));

				double mapRef_x = rotate_x + ref_x;
				double mapRef_y = rotate_y + ref_y;

				next_x_vals.push_back(mapRef_x);
				next_y_vals.push_back(mapRef_y);

			}

		



          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
#ifdef UWS_VCPKG
			// code fixed for latest uWebSockets
			ws->send(msg.data(), msg.length(), uWS::OpCode::TEXT);
#else
			// leave original code here
			ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
#endif
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
#ifdef UWS_VCPKG
		// code fixed for latest uWebSockets
		ws->send(msg.data(), msg.length(), uWS::OpCode::TEXT);
#else
		// leave original code here
		ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
#endif
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

#ifdef UWS_VCPKG
  // code fixed for latest uWebSockets
  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> *ws, uWS::HttpRequest req) {
#else
  // leave original code here
  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
#endif
	  std::cout << "Connected!!!" << std::endl;
  });

#ifdef UWS_VCPKG
  // code fixed for latest uWebSockets
  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> *ws, int code, char *message, size_t length) {
#else
  // leave original code here
  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
#endif
#ifdef UWS_VCPKG
	  // code fixed for latest uWebSockets
	  ws->close();
#else
	  // leave original code here
	  ws.close();
#endif
	  std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen("127.0.0.1", port)){
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
