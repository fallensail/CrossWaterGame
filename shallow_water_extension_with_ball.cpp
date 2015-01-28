/**
 * @file shallow_water_extension.cpp
 * Implementation of a shallow water system using Mesh
 * Added vertical force for floating obejcts forces
 * Added source functions for uneven water surface
 * Combined with Graph for game ending visualization
 * ---- Combined with Meshed mass spring for game winning visualization
 * 
 * @brief Reads in two files specified on the command line and two integer
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 * Integer: Initial conditions. 0 for static
 * Example: ./shallow_water data/tub3.* 0
 */

#include <fstream>
#include <cmath>
#include <string>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <queue>
#include <ctime>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Point.hpp"
#include "Mesh.hpp"
#include "Graph.hpp"

#include "CollisionDetector.hpp"


// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

// This is one of the example object that can be floating on the water
// It should contain weights, length, its location and speed available
// Again this is just a simple example of floating objects which can not be shown in viewer
// Main fuctionality is generalized using template in the latter part for genaric usage
// Small ship object used in the game provided
struct Ship {
  double Weight;
  double length;
  double center_x;
  double center_y;
  double speed_x;
  double speed_y;
  double Area;

  /** Default constructor for ship
   *  Create a generic ship that has weight, length, center_x, center_y, speed_x, speed_y set as default value
   *  Area can be changed and is calculated depending on the user need
   **/
  Ship(): Weight(0.015), length(0.15), center_x(-0.8), center_y(-0.9), speed_x(0), speed_y(0) {
    double e = length/2;
    Area = 3.1416 * e * e;
  }

  /** Constructor for ship
   * @param[in] double Weight
   * @param[in] double length
   * @param[in] double center_x
   * @param[in] double center_y
   * @param[in] double speed_x
   * @param[in] double speed_y
   *
   * Create a generic ship that has weight, length, center_x, center_y, speed_x, speed_y
   * Area can be changed and is calculated depending on the user need
   * 
   **/
  Ship(double w, double a, double c_x, double c_y, double s_x, double s_y)
      : Weight(w), length(a), center_x(c_x), center_y(c_y), speed_x(s_x), speed_y(s_y) {
    double e = length/2;
    Area = 3.1416 * e * e;
  }

  /** Check if a point is under the ship area
   * @param[in] Point a point position to check
   * @param[out] bool true if it is under the ship, false if it is not
   * 
   * This function check if a point is under the ship. 
   * This can be changed based on the shape of the ship and direction of the ship
   */
  bool checkCoverNode(Point p) {
    if (norm(p - Point(center_x, center_y, 0))
       < sqrt(Area/3.1416)) {
       return true;
    } else {
      return false;
    }
  }

  /** Move the ship location in time dt
   * @param[in] double dt small change in time period
   * 
   * This fucntion changes the ship location in short time dt based on 
   * current speed current location
   */
  void move(double dt) {
    center_x += speed_x * dt;
    center_y += speed_y * dt;
  }
};

 /** Source calculcator of partial derivative of x by varying position x and y
   * @param[in] double x
   * @param[in] double y
   * @param[out] double partial derivative of x at position x and y
   *
   * This function can be changed to tailer different shapes of source
   */
double dx_value(double x, double y){
  if (x > 0) {
    x = -x;
  }
  return 0.12*2*x+0*y;
};

/** Source calculcator of partial derivative of y by varying position x and y
   * @param[in] double x
   * @param[in] double y
   * @param[out] double partial derivative of y at position x and y
   *
   * This function can be changed to tailer different shapes of source
   */
double dy_value(double x, double y){
  if (y > 0) {
    y = -y;
  }
  return 0.15*2*y+0*x;
};

/** @struct Mesh::NodeData
  * information associated with node
  */    
struct NodeData {
  QVar Q;
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
  NodeData():Q(0.0,0.0,0.0){}
};

/** @struct Mesh::TriData
  * information associated with triangles
  */
struct TriData {
  QVar Q;  //Qk: the average value of Q inside this triangle
  std::vector<QVar> F; //The transition over an edge
  Point n; //the outward surface normal vector
  TriData():Q(QVar()), F(3, QVar()){}
};

/** @struct Mesh::EdgeData
  * information associated with edge
  */
struct EdgeData {
  double L;
  double K;
};

typedef Mesh<NodeData, EdgeData, TriData> MeshType;
typedef typename MeshType::Triangle Triangle;
typedef typename MeshType::node_type Node;
typedef typename MeshType::edge_type Edge;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the {n+1} node positions
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    n.position() += n.value().velocity * dt;
  }
  // Compute the {n+1} node velocities
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    n.value().velocity += force(n, t) * (dt / n.value().mass);
  }
  return t + dt;
}

//Force Function to calculate gravity
struct GravityForce
{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return (Point(0, 0, -grav)*n.value().mass);
  }
};

//Force Function to calculate spring force
struct MassSpringForce
{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point spring = Point(0, 0, 0);  //spring force
    //add up all spring forces
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      auto incident = *it;

      if(incident.node1()==n)
        spring += ((incident.value().K)*(incident.node2().position()-incident.node1().position())/incident.length()*(incident.length()-incident.value().L));
      else
        spring += ((incident.value().K)*(incident.node1().position()-incident.node2().position())/incident.length()*(incident.length()-incident.value().L));
    }
    (void) t;
    return spring;
  }
};

//Force Function to calculate damp force
struct DampingForce
{
  DampingForce(double coef): c(coef) {}
  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    return (-(c*n.value().velocity));
  }
  double c;
};


// The wind force
struct WindForce {
  WindForce(Point wind): w(wind) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    double c = 0.00004;
    auto normal = Point(0,0,0);
    for (auto it=n.triangle_begin(); it!=n.triangle_end(); ++it){
      Point tnorm;
      tnorm = cross((*it).node(0).position()-(*it).node(1).position(), 
        (*it).node(0).position()-(*it).node(2).position());
      tnorm = tnorm/norm(tnorm);
      normal = normal + tnorm;
    }
    (void) t;
    return c*dot((w-n.value().velocity),normal)*normal;
  }
  Point w;
};


//the air pressure force
template <typename NODE, typename G>
struct PressureForce
{
  PressureForce(double p_out, double c, G* graph): P_out(p_out), C(c), g(graph) {}

  Point operator()(NODE n, double t) {

    //if n.index()==0, update the volume, center and normal vector,
    //P_diff, for each node
    if(n.index() == 0){
      //update the center
      center = Point(0, 0, 0);
      for (auto it=g->node_begin(); it != g->node_end(); ++it){
        center += (*it).position()/g->num_nodes();
      }

      //update the outward normal vector
      for (auto it=g->triangle_begin(); it != g->triangle_end(); ++it){
        Point tnorm;
        tnorm = cross((*it).node(0).position()-(*it).node(1).position(), 
          (*it).node(0).position()-(*it).node(2).position());
        if (dot(tnorm, center-(*it).node(0).position())>0)
          tnorm = -tnorm;
        tnorm = tnorm/norm(tnorm);
        (*it).value().n = tnorm;
      }
      //update the Volume
      V = 0;
      for (auto it=g->triangle_begin(); it != g->triangle_end(); ++it){
        V += (*it).value().n.z*(*it).area()*((*it).node(0).position().z +
          (*it).node(1).position().z + (*it).node(2).position().z)/3;
      }

      //P_diff
      P_diff = C/V - P_out;
    }

    //for any node, calculate the force
    Point force = Point(0, 0, 0);
    for (auto it=n.triangle_begin(); it != n.triangle_end(); ++it){
      force += P_diff*(*it).area()*(*it).value().n/3;
    }
    (void) t;
    //std::cout<<force<<std::endl;
    return force;
  }
private:
  double P_out; //the pressure outside the ball
  double C;  //nRT
  double V;  //the volumn of the ball
  double P_diff; //P_inside-P_out
  Point center; //The point in the center of the ball
  G* g;
};

//Force function which represents the combined effects of F1 and F2
template<typename F1,typename F2>
struct CombinedForce{
  F1 force1;
  F2 force2;
  CombinedForce(F1 f1=F1(), F2 f2=F2()):force1(f1), force2(f2){}
  template <typename NODE>
  Point operator() (NODE n, double t){
    return (force1(n, t)+force2(n, t));
  }
};

//Combine the effects of two forces
template<typename F1,typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1 = F1(), F2 f2 = F2()){
  return CombinedForce<F1, F2>(f1, f2);
}
//Combine the effects of three forces
template<typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1, F2>, F3> make_combined_force(F1 force1, F2 force2, F3 force3){
  return make_combined_force(make_combined_force(force1, force2), force3);
}


/** Function object for calculating shallow-water flux.
 *          |n
 *   T_k    |---> n = (nx,ny)   T_m
 *   QBar_k |                   QBar_m
 *          |
 * @param[in] nx, ny Defines the 2D outward normal vector n = (@a nx, @a ny)
 *            from triangle T_k to triangle T_m. The length of n is equal to the
 *            the length of the edge, |n| = |e|.
 * @param[in] dt The time step taken by the simulation. Used to compute the
 *               Lax-Wendroff dissipation term.
 * @param[in] qk The values of the conserved variables on the left of the edge.
 * @param[in] qm The values of the conserved variables on the right of the edge.
 * @return The flux of the conserved values across the edge e
 */
template <typename OBJ>
struct EdgeFluxCalculator {
  QVar operator()(double nx, double ny, double dt,
                  const QVar& qk, const QVar& qm, Point p1, Point p2, std::vector<OBJ>& obj_vector) {
    // Normalize the (nx,ny) vector
    double n_length = std::sqrt(nx*nx + ny*ny);
    nx /= n_length;
    ny /= n_length;

    // The velocities normal to the edge
    double wm = (qm.hx*nx + qm.hy*ny) / qm.h;
    double wk = (qk.hx*nx + qk.hy*ny) / qk.h;

    // Lax-Wendroff local dissipation coefficient
    double vm = sqrt(grav*qm.h) + sqrt(qm.hx*qm.hx + qm.hy*qm.hy) / qm.h;
    double vk = sqrt(grav*qk.h) + sqrt(qk.hx*qk.hx + qk.hy*qk.hy) / qk.h;
    double a  = dt * std::max(vm*vm, vk*vk);

    // Helper values
    double scale = 0.5 * n_length;
    double gh2   = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h);

    for (unsigned i = 0; i < obj_vector.size(); ++i) {
      if (obj_vector[i].checkCoverNode(p1)) {
        gh2 = gh2  +  obj_vector[i].Weight*qm.h/1/obj_vector[i].Area;
      }
      if (obj_vector[i].checkCoverNode(p2)) {
        gh2 = gh2 + obj_vector[i].Weight*qk.h/1/obj_vector[i].Area;
      }
    }
    // Simple flux with dissipation for stability
    return QVar(scale * (wm*qm.h  + wk*qk.h)           - a * (qm.h  - qk.h),
                scale * (wm*qm.hx + wk*qk.hx + gh2*nx) - a * (qm.hx - qk.hx),
                scale * (wm*qm.hy + wk*qk.hy + gh2*ny) - a * (qm.hy - qk.hy));
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) {
    // Change this to plot something other than the
    // positions of the nodes
    return Point(n.position().x, n.position().y, n.value().Q.h);
  }
};


/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX, typename OBJ>
double hyperbolic_step(MESH& m, FLUX& f, double t, double dt, std::vector<OBJ>& obj_vector) {
  // Step the finite volume model in time by dt.
  // Pseudocode:
  // Compute all fluxes. (before updating any triangle Q_bars)
  // For each triangle, update Q_bar using the fluxes as in Equation 8.
  //  NOTE: Much like symp_euler_step, this may require TWO for-loops
  for (auto i = m.edge_begin(); i != m.edge_end(); ++i) {

    if ((*i).triangle1().index() != (unsigned) -1 && (*i).triangle2().index() != (unsigned) -1) {
      MeshType::Triangle trik = (*i).triangle1();
      MeshType::Triangle trim = (*i).triangle2();
      unsigned int edge_k = 0;
      unsigned int edge_m = 0;
      //which edge (*i) is in trik and trim
      while(trik.node(edge_k).index()== (*i).node1().index() 
        || trik.node(edge_k).index()== (*i).node2().index() )
        ++edge_k;
      while(trim.node(edge_m).index()== (*i).node1().index() 
        || trim.node(edge_m).index()== (*i).node2().index() )
        ++edge_m;
      QVar flux = f(trik.normal(edge_k).x, trik.normal(edge_k).y, dt, trik.value().Q, trim.value().Q, (*i).node1().position(), (*i).node2().position(), obj_vector);
      trik.value().F[edge_k] = flux;
      trim.value().F[edge_m] = -flux;
    } else {
      MeshType::Triangle trik;
      if ((*i).triangle1().index() != (unsigned) -1)
        trik = (*i).triangle1();
      else
        trik = (*i).triangle2();
      unsigned int edge_k = 0;
      while(trik.node(edge_k).index()== (*i).node1().index() 
        || trik.node(edge_k).index()== (*i).node2().index() )
        ++edge_k;
      QVar flux = f(trik.normal(edge_k).x, trik.normal(edge_k).y, dt, trik.value().Q, QVar(trik.value().Q.h, 0, 0), (*i).node1().position(), (*i).node2().position(), obj_vector);
      trik.value().F[edge_k] = flux;
    }
  }

  for(auto i = m.triangle_begin(); i != m.triangle_end(); ++i){
    QVar sum = QVar(0, 0, 0);
    Point center = Point(0,0,0);
    for (int j = 0; j < 3; ++j){
      sum += (*i).value().F[j];
      center += (*i).node(j).position();
    }
    center = center/3;
    (*i).value().Q = (*i).value().Q-dt/(*i).area()*sum + 
      dt * QVar(0, -grav*(*i).value().Q.h*dx_value(center.x,center.y), 
                   -grav*(*i).value().Q.h*dy_value(center.x,center.y));
  }
  return t + dt;
}

/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
  // Translate the triangle-averaged values to node-averaged values
  // Implement Equation 9 from your pseudocode here
  for (auto it = m.node_begin(); it != m.node_end(); ++it){
    double sumarea=0;
    QVar sumQ = QVar(0, 0, 0);
    for(auto j = m.triangle_begin(*it); j != m.triangle_end(*it); ++j){
      sumarea += (*j).area();
      sumQ += (*j).value().Q * (*j).area();
    }
    (*it).value().Q = sumQ/sumarea;
    //(*it).position().z = (*it).value().Q.h;
  }
}

// Construct a Color functor for ship and water and view with the SDLViewer
template <typename Ship>
struct ColorFunctor {
	std::vector<Ship>& obj_vector;
	
  template <typename NODE>
	CS207::Color operator()(const NODE& n) const {
	  if (norm(n.position() - Point(0.75,1,0)) < 0.1) {
	    return CS207::Color(0.9,0.2,0.8);
	  }
    bool check = false;
		if (obj_vector.size() >= 1 && obj_vector[0].checkCoverNode(n.position())) {
			return CS207::Color(0.1,0.9,0.1);
		}
		if (obj_vector.size() >= 2) {
      for (unsigned i = 1; i < obj_vector.size(); ++i) {
        if (obj_vector[i].checkCoverNode(n.position())) {
          check = true;
        }
      }
		}
    if (check) {
      return CS207::Color(0.9,0.1,0.2);
    }
		double h = n.value().Q.h;
		if (h < 1.01) {
			return CS207::Color(0.0,0.0,1.0);
		} else if (h >= 1.2) {
			return CS207::Color(1.0,1.0,1.0);
		} else {
			return CS207::Color((h-1)*4.8, (h-1)*4.8, 1.0);
		}
	};
};

// Construct a Color functor for water only and view with the SDLViewer
template <typename Ship>
struct ColorFunctorWater {
  std::vector<Ship>& obj_vector;
  
  template <typename NODE>
  CS207::Color operator()(const NODE& n) const {
    double h = n.value().Q.h;
    if (h < 1.01) {
      return CS207::Color(0.0,0.0,1.0);
    } else if (h >= 1.2) {
      return CS207::Color(1.0,1.0,1.0);
    } else {
      return CS207::Color((h-1)*4.8, (h-1)*4.8, 1.0);
    }
  };
};

// Construct a Color functor for sphere ball only and view with the SDLViewer
struct ColorFunctorBall {

  template <typename NODE>
  CS207::Color operator()(const NODE& n) const {
    double h = n.position().z;
    if (h < 0.9) {
      return CS207::Color(0.0,0.0,1.0);
    } else if (h >= 1.1) {
      return CS207::Color(1.0,1.0,1.0);
    } else {
      return CS207::Color((h-0.9)*5, (h-0.9)*5, 1.0);
    }
  };
};


// Construct a Color functor for graph skull and view with the SDLViewer
struct ColorFn {
  int max_;
  ColorFn(int max) {
    max_ = max;
  }
  
  CS207::Color operator() (Graph<int>::Node n) {
    float fraction = 1.0-((float)n.value())/(max_+1);
    if (fraction > 1) {
      fraction = 1;
    } else if (fraction < 0) {
      fraction = 0;
    }
    return CS207::Color::make_heat(fraction);
  }
};

// Listener for user controled ship that can be moved based on user arrow key input
struct Listener_MyShip: public CS207::SDLViewer::Listener{
  Ship& ship;
  Listener_MyShip(Ship& s): ship(s){}
  void handle(SDL_Event e){
    switch (e.type) {

      case SDL_KEYDOWN: {
        // Keyboard 'arrow right' to increase wind
        if (e.key.keysym.sym == SDLK_RIGHT){
          ship.center_y += 0.05;
        }
        // Keyboard 'arrow left' to decrease wind
        if (e.key.keysym.sym == SDLK_LEFT){
          ship.center_y -= 0.05;
        }
        // Keyboard 'arrow up' to increase location of wind
        if (e.key.keysym.sym == SDLK_UP){
          ship.center_x -= 0.05;
        }
        // Keyboard 'arrow down' to decrease location of wind
        if (e.key.keysym.sym == SDLK_DOWN){
          ship.center_x += 0.05;
        }

        // Constrain the user boat to be within the boundary
        if (ship.center_x > 1) ship.center_x = 1;
        if (ship.center_x < -1) ship.center_x = -1;
        if (ship.center_y > 1) ship.center_y = 1;
        if (ship.center_y < -1) ship.center_y = -1;
      } break;
    }
  }
};

/** Comparator that compares the distance from a given point p.
 */
struct MyComparator {
   Point p_;
   MyComparator(const Point& p) : p_(p) {
   };

   template <typename NODE>
   bool operator()(const NODE& node1, const NODE& node2) const {
    Point::size_type dist1 = (node1.position().x - p_.x)*(node1.position().x - p_.x) 
      + (node1.position().y - p_.y)*(node1.position().y - p_.y)
      + (node1.position().z - p_.z)*(node1.position().z - p_.z);
    Point::size_type dist2 = (node2.position().x - p_.x)*(node2.position().x - p_.x) 
      + (node2.position().y - p_.y)*(node2.position().y - p_.y)
      + (node2.position().z - p_.z)*(node2.position().z - p_.z);
    return dist1<dist2;
  }
};


/** Calculate shortest path lengths in @a g from the nearest node to @a point.
 * @param[in,out] g Input graph
 * @param[in] point Point to find the nearest node to.
 * @post Graph has modified node values indicating the minimum path length
 *           to the nearest node to @a point
 * @post Graph nodes that are unreachable to the nearest node to @a point have
 *           the value() -1.
 * @return The maximum path length found.
 *
 * Finds the nearest node to @a point and treats that as the root node for a
 * breadth first search.
 * This sets node's value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(Graph<int>& g, const Point& point) {
  MyComparator pred = MyComparator(point);
  auto root_node = *std::min_element(g.node_begin(), g.node_end(), pred);
  
  int max = 0;
  
  for (Graph<int>::NodeIterator iter = g.node_begin(); iter != g.node_end(); ++iter) {
    (*iter).value() = -1;
  }
  root_node.value() = 0;
  
  // BFS:
    std::queue<Point::size_type> Q;
 
    /** Keeps track of explored vertices */
    std::vector<bool> explored;
 
    /** Initailized all vertices as unexplored */
    for (Point::size_type i = 0; i < g.num_nodes(); ++i)
      explored.push_back(false);
 
    /** Push initial vertex to the queue */
    Q.push(root_node.index());
    explored[root_node.index()] = true; /** mark it as explored */

    /** Unless the queue is empty */
    while (!(Q.size() == 0)) {
        /** Pop the vertex from the queue */
        Point::size_type v = Q.front();
        Q.pop();

        /** From the explored vertex v try to explore all the
        connected vertices */
        for (Graph<int>::IncidentIterator iter = g.node(v).edge_begin(); g.node(v).edge_end() != iter; ++iter) {
            /** Explores the vertex if it is connected to v
            and if it is unexplored */
            if (!explored[(*iter).node2().index()]) {
                /** adds the new vertex to the queue */
                Q.push((*iter).node2().index());
                (*iter).node2().value() = (*iter).node1().value() + 1;
                if (max < (*iter).node2().value()) {
                  max = (*iter).node2().value();
                }
                /** marks the new vertex as visited */
                explored[(*iter).node2().index()] = true;
            }
        }
    }
  return max;
}


int main(int argc, char* argv[])
{
  srand(time(NULL));
  // Check arguments
  if (argc < 4) {
    std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE NUM_CONDITION\n";
    exit(1);
  }

  // Construct a Graph for game over scene
  typedef Graph<int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file2("data/large.nodes");
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p2;
  while (CS207::getline_parsed(nodes_file2, p2))
    nodes.push_back(graph.add_node(p2));

  // Create a tets_file from the second input argument
  std::ifstream tets_file2("data/large.tets");
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t2;
  while (CS207::getline_parsed(tets_file2, t2))
    for (unsigned i = 1; i < t2.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t2[i]], nodes[t2[j]]);

  // Use shortest_path_lengths to set the node values to the path lengths
  int max = shortest_path_lengths(graph, Point(-1,0,1));

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Construct a Mesh sphere ball for game winning scene
  MeshType mesh2;
  std::vector<typename MeshType::node_type> mesh_node2;
  // Read all Points and add them to the Mesh
  std::ifstream nodes_file3("data/sphere200.nodes");
  Point p3;
  while (CS207::getline_parsed(nodes_file3, p3)) {
    mesh_node2.push_back(mesh2.add_node(p3));
  }

  // Read all mesh triangles and add them to the Mesh sphere ball
  std::ifstream tris_file3("data/sphere200.tris");
  std::array<int,3> t3;
  while (CS207::getline_parsed(tris_file3, t3)) {
    mesh2.add_triangle(mesh_node2[t3[0]], mesh_node2[t3[1]], mesh_node2[t3[2]]);
  }

  std::cout << mesh2.num_nodes() << " "
            << mesh2.num_edges() << " "
            << mesh2.num_triangles() << std::endl;

  // initialize values for Mesh sphere ball
  for (auto it = mesh2.node_begin(); it != mesh2.node_end(); ++it){
    (*it).value().mass = float(1)/mesh2.num_nodes();
    (*it).value().velocity = Point(0, 0, 0);
  }

  for (auto it = mesh2.node_begin(); it != mesh2.node_end(); ++it)
  {
    for (auto j = (*it).edge_begin(); j != (*it).edge_end(); ++j){
       (*j).value().L = (*j).length();
       (*j).value().K = 400;
    }
  }

  for (auto it = mesh2.node_begin(); it != mesh2.node_end(); ++it){
    (*it).value().mass = float(1)/mesh2.num_nodes();
    (*it).value().velocity = Point(0, 0, 0);
  }

  for (auto it = mesh2.node_begin(); it != mesh2.node_end(); ++it){
    (*it).position().z += 3;
  }
  
  // create Mesh shallow water
  MeshType mesh;
  std::vector<typename MeshType::node_type> mesh_node;

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
      mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;

  // Set the initial conditions of mesh shallow water
  if (strcmp(argv[3], "1") == 0) {
    for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) {
      (*it).value().Q = QVar(1-0.75*exp(-80*((pow((*it).position().x-0.75, 2.0)+(*it).position().y*(*it).position().y))), 0, 0);
    }
  } else if (strcmp(argv[3], "2") == 0) {
    for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) {
      if (pow(((*it).position().x-0.75),2) + (*it).position().y*(*it).position().y -0.15*0.15< 0) { 
        (*it).value().Q = QVar(1.75,0,0);
      } else {
        (*it).value().Q = QVar(1,0,0);
      }
    }
  } else if (strcmp(argv[3], "3") == 0) {
    for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) {
      if ((*it).position().x < 0) { 
        (*it).value().Q = QVar(1.75,0,0);
      } else {
        (*it).value().Q = QVar(1,0,0);
      }
    }
  } else if (strcmp(argv[3],"0") == 0) {
    for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) {
        (*it).value().Q = QVar(1,0,0);
    }
    std::cout << "No initial condition for sepcifying 0 but OK." << std::endl;
  } else {
    std::cerr << "SEPCIFY INITIAL CONDITION BY ADDING THE FOURTH VALUE INT 0-2\n";
    exit(1);
  }
  
  for (auto it=mesh.triangle_begin(); it!=mesh.triangle_end(); ++it) {
    (*it).value().Q = ((*it).node(0).value().Q + (*it).node(1).value().Q + (*it).node(2).value().Q)/3;
  }

  // Launch the SDLViewer
  CS207::SDLViewer viewer;

  // add water into viewer
  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

  viewer.center_view();
  viewer.launch();

  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt
  double min_edge_length = (*mesh.edge_begin()).length();
  
  for (auto iter = mesh.edge_begin(); iter != mesh.edge_end(); ++iter) {
    if (min_edge_length > (*iter).length()) {
      min_edge_length = (*iter).length();
    }
  }

  double max_height = mesh.node(0).value().Q.h;

  for (auto iter = mesh.node_begin(); iter != mesh.node_end(); ++iter) {
    if (max_height < (*iter).value().Q.h) {

      max_height = (*iter).value().Q.h;
    }
  }
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_height));
  double t_start = 0;
  double t_end = 10;

  // Preconstruct a Flux functor
  EdgeFluxCalculator<Ship> f;
  
  Ship obj0 = Ship();
  
  std::cout << dt << std::endl;

  std::vector<Ship> obj_vector;
  obj_vector.push_back(obj0);
  for (unsigned j = 0; j < 6; ++j) {
    obj_vector.push_back(Ship(0.03, 0.2, 1 - 2*(rand()%100/100.0), -0.5 + j*0.25, -2.0 - 1.5*(rand()%100/100.0), 0));
  }
  
  ColorFunctor<Ship> ColFn {obj_vector};
  ColorFunctorWater<Ship> ColFnWater {obj_vector};
  ColorFunctorBall ColFnBall;
  ColorFn cf = ColorFn(max);

  // new for viewer user ship interaction of moving user ship
  std::shared_ptr<CS207::SDLViewer::Listener> ptr_ship(new Listener_MyShip(obj_vector[0]));
  viewer.add_listener(ptr_ship);

  // force creation for Meshed ball force
  WindForce wind_force(Point(0,0,0));
  PressureForce<typename MeshType::node_type, MeshType> pressure_force(0.2, 1, &mesh2);
  DampingForce damp_force(float(1)/mesh2.num_nodes());
  auto force = make_combined_force(MassSpringForce(), GravityForce(), make_combined_force(pressure_force, damp_force, wind_force));

  bool winning = false; // winning check

  // Start time-step
  for (double t = t_start; t < t_end; t += dt) {

    // random fire ship speed and presenting time generation
    for (unsigned i = 1; i < obj_vector.size(); ++i) {
      obj_vector[i].move(dt);
      obj_vector[i].center_x = fmod(obj_vector[i].center_x - 1, 2.0 + 1.5*(rand()%100/100.0)) + 1;
      obj_vector[i].speed_x += dt*(rand()%100/100.0 - 0.3);
    }

    // dead condtion check
    bool checker = false;
    for (unsigned i = 1; i < obj_vector.size(); ++i) {
      if (norm(Point(obj_vector[0].center_x, obj_vector[0].center_y, 0) - Point(obj_vector[i].center_x, obj_vector[i].center_y, 0)) < 
        obj_vector[0].length/2 + obj_vector[i].length/2) {
        checker = true;
      }
    }

    // Lost Scene
    if (checker || t > (9.9 - dt)) {
      auto node_map = viewer.empty_node_map(graph);

      for (auto iter = graph.node_begin(); iter != graph.node_end(); ++iter) {
        (*iter).position() += Point(0,0,1.5);
        (*iter).position() = Point(-(*iter).position().y, (*iter).position().x, (*iter).position().z);
      }

      viewer.add_nodes(graph.node_begin(), graph.node_end(), cf, node_map);
    //MyNodePredicate predicate;
    //viewer.add_nodes(make_filtered(graph.node_begin(), graph.node_end(), predicate), make_filtered(graph.node_end(), graph.node_end(), predicate), cf, node_map);
      viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
      return 0;
    }

    // Winning Scene Start
    if (norm(Point(obj_vector[0].center_x, obj_vector[0].center_y, 0) - Point(0.75,1,0)) < 0.1 && !winning) {
    //if (t > 10*dt && !winning) {
      winning = true;
      obj_vector.clear();
      obj_vector.push_back(Ship(0.05, 0.5, 0.0, 0.0, 0.0, 0.0));

      viewer.add_nodes(mesh2.node_begin(), mesh2.node_end(), ColFnWater, node_map);
      viewer.add_edges(mesh2.edge_begin(), mesh2.edge_end(), node_map);
    }

    // Step forward in time with forward Euler for meshed shallow water
    hyperbolic_step(mesh, f, t, dt, obj_vector);
    // Update node values with triangle-averaged values
    post_process(mesh);

    // Update the viewer with new node positions depending on winning already or not yet winning
    if (!winning) {
      viewer.add_nodes(mesh.node_begin(), mesh.node_end(), ColFn, NodePosition(), node_map);
    } else {
      viewer.add_nodes(mesh.node_begin(), mesh.node_end(), ColFnWater, NodePosition(), node_map);
    }
    viewer.set_label(t);

    // winning scene simulation
    if (winning) {
      obj_vector[0] = Ship(0, 0.5, 0.0, 0.0, 0.0, 0.0);

      // add Meshed sphere ball
      viewer.add_nodes(mesh2.node_begin(), mesh2.node_end(), ColFnBall, node_map);

      // set shallow water z axis poision ready for collision detection
      for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it){
        (*it).position().z = (*it).value().Q.h;
      }
      
      // moving Meshed sphere ball
      symp_euler_step(mesh2, t, dt, force);

      // collision detection check
      CollisionDetector<MeshType> c;
      c.add_object(mesh);
      c.add_object(mesh2);
      c.check_collisions();
      std::vector<unsigned> collision;
      std::vector<unsigned> collision2;

      // store collision info
      unsigned collision_counter = 0;
      for (auto it=c.begin(); it!= c.end(); ++it){
        auto boom = *it;
        if (boom.mesh1 == &mesh)
          ++collision_counter;
      }

      if (collision_counter > 3) {
        obj_vector[0] = Ship(1.2, 0.5 + 0.01 * collision_counter, 0.0, 0.0, 0.0, 0.0);
        WindForce wind_force2(Point(0,0, collision_counter));
        force = make_combined_force(MassSpringForce(), 
          GravityForce(), make_combined_force(pressure_force, damp_force, wind_force2));
      }

      for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it){
        (*it).position().z = 0;
      }
      std::cout<<"Collision Detection (number of nodes): " << collision_counter << std::endl;
    }

    // These lines slow down the animation for small meshes.
    // Feel free to remove them or tweak the constants.
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.05);
  }
  return 0;
}
