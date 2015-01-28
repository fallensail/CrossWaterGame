/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;
// damping constant
double c;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
};

// HW2 #1 YOUR CODE HERE
// Define your Graph type
typedef Graph<NodeData, double> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;

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

  //std::cout<<"symp_euler_step started (inside)"<<std::endl;
  // Compute the {n+1} node positions
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().velocity * dt;
  }
  
  // Compute the {n+1} node velocities
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().velocity += force(n, t) * (dt / n.value().mass);
  }
  
  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force being applied to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  Point operator()(Node n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;
    
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    }
    
    Point force = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto a = n.position() - (*it).node2().position();
      double b = (*it).length();
      double c = (*it).length()- (*it).value();
      force = force + (-100)*(a/b)*c;
    }
    
    force = force + n.value().mass*Point(0,0, -grav);
    return force;
  }
};

// Force Functors 
/** Gravity force applied on an object mass*Point(0,0,-grav)*/
struct GravityForce {
  /** Calculate gravity force
   * @param[in]     n      Node
   * @param[in]     t      double time
   * @return a Point object as a force
   *
   */
  Point operator() (Node n, double t) {
    (void) t;
    return n.value().mass*Point(0,0, -grav);
  }
};

/** spring force applied on an object
 */
struct MassSpringForce {
  /** Calculate spring force
   * @param[in]     n      Node
   * @param[in]     t      double time
   * @return a Point object as a force
   */
  Point operator() (Node n, double t) {
    (void) t;
    Point force = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto a = n.position() - (*it).node2().position();
      double b = (*it).length();
      double c = (*it).length()- (*it).value();
      force = force + (-100)*(a/b)*c;
    }
    return force;
  }
};

/** damping force applied on an object
 */
struct DampingForce {
  /** Calculate spring force
   * @param[in]     n      Node
   * @param[in]     t      double time
   * @return a Point object as a force
   */
  Point operator() (Node n, double t) {
    (void) t;
    Point force = Point(0,0,0);
    return force - c*n.value().velocity;
  }
};

/** MetaForce to combine/aggregate any two force
 *
 * @tparam F1 force functor 1
 * @tparam F2 force functor 2
 */
template<typename F1, typename F2>
class MetaForce {
  public:
    F1 f1;
    F2 f2;
    /** MetaForce to combine/aggregate any two force
     *
     * @param[in]     n      Node
     * @param[in]     t      double time
        * @return a Point object as a force
     */
    Point operator() (Node n, double t) {
      return f1(n, t) + f2(n, t);
    }
};

/** Create a new MetaForce that combines two other forces
 * @param[in]     f1  Force functor 1
 * @param[in]     f2  Force functor 2
 * @return a new MetaForce
 *
 * @tparam F1 force functor 1
 * @tparam F2 force functor 2
 */
template<typename F1, typename F2>
MetaForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
  MetaForce<F1, F2> mf;
  mf.f1 = f1;
  mf.f2 = f2;
  return mf;
}

/** Create a new MetaForce that combines three other forces
 * @param[in]     f1  Force functor 1
 * @param[in]     f2  Force functor 2
 * @param[in]     f2  Force functor 3
 * @return a new MetaForce
 *
 * @tparam F1 force functor 1
 * @tparam F2 force functor 2
 * @tparam F2 force functor 3
 */
template<typename F1, typename F2, typename F3>
MetaForce<MetaForce<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  MetaForce<F1, F2> mf;
  mf.f1 = f1;
  mf.f2 = f2;
  MetaForce<MetaForce<F1, F2>, F3> mf2;
  mf2.f1 = mf;
  mf2.f2 = f3;
  return mf2;
}

// Constraint functors
/** Specific constraint applied on specific node*/
struct SpecificConstraint {
  /** Change node position and velocity based on constraint
   * @param[in]     g      graph
   * @param[in]     t      double time
   */
  void operator() (GraphType& g, double t){
    (void) t;
    for (auto iter = g.node_begin(); iter !=g.node_end(); ++iter) {
      if ((*iter).position() == Point(0,0,0) || (*iter).position() == Point(1,0,0)) {
        (*iter).value().velocity = Point(0,0,0);
      }
    }
  }
};

// Constraint functors
/** Plane constraint applied on specific node*/
struct PlaneConstraint {
  double plane_z = -0.75;
  /** Change node position and velocity based on constraint
   * @param[in]     g      graph
   * @param[in]     t      double time
   */
  void operator() (GraphType& g, double t) {
    (void) t;
    for (auto iter = g.node_begin(); iter !=g.node_end(); ++iter) {
      if ((*iter).position().z < plane_z) {
        (*iter).position().z = plane_z;
        (*iter).value().velocity.z = 0;  
      }
    }
  }
};

// Constraint functors
/** Sphere constraint applied on specific node*/
struct SphereConstraint {
  double radius = 0.15;
  Point center = Point(0.5, 0.5, -0.5);
  
  /** Change node position and velocity based on constraint
   * @param[in]     g      graph
   * @param[in]     t      double time
   */
  void operator() (GraphType& g, double t) {
    (void) t;
    for (auto iter = g.node_begin(); iter !=g.node_end(); ++iter) {
      double dist = norm((*iter).position() - center);
      if (dist < radius) {
        Point norm = ((*iter).position() - center)/dist;
        
        (*iter).position() = norm*radius + center;
        (*iter).value().velocity = (*iter).value().velocity - dot((*iter).value().velocity, norm)*norm; 
      }
    }
  }
};

// Constraint functors
/** Sphere constraint applied on specific node*/
struct RemoveSphereConstraint {
  double radius = 0.15;
  Point center = Point(0.5, 0.5, -0.5);
  
  /** Change node position and velocity based on constraint
   * @param[in]     g      graph
   * @param[in]     t      double time
   */  
  void operator() (GraphType& g, double t) {
    (void) t;
    for (auto iter = g.node_begin(); iter !=g.node_end(); ++iter) {
      double dist = norm((*iter).position() - center);
      if (dist < radius) {
        g.remove_node((*iter));
      }
    }
  }  
};

/** MetaConstraint to combine/aggregate any two constraint
 *
 * @tparam C1 force constraint 1
 * @tparam C2 force constraint 2
 */
template<typename C1, typename C2>
class MetaConstraint {
  public:
    C1 c1;
    C2 c2;
    /** MetaForce to combine/aggregate any two constraint
     *
     * @param[in]     n      Node
     * @param[in]     t      double time
     * @return a Point object as a constraint
     */
    void operator() (GraphType& g, double t) {
      c1(g, t);
      c2(g, t);
    }
};

/** Create a new MetaConstraint that combines two other constraints
 * @param[in]     c1  Constraint  functor 1
 * @param[in]     c2  Constraint  functor 2
 * @return a new MetaForce
 *
 * @tparam C1 Constraint  functor 1
 * @tparam C2 Constraint  functor 2
 */
template<typename C1, typename C2>
MetaConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
  MetaConstraint<C1, C2> mc;
  mc.c1 = c1;
  mc.c2 = c2;
  return mc;
}

/** Create a new MetaConstraint that combines three other constraints
 * @param[in]     c1  Constraint  functor 1
 * @param[in]     c2  Constraint  functor 2
 * @param[in]     c3  Constraint  functor 3
 * @return a new MetaForce
 *
 * @tparam C1 Constraint  functor 1
 * @tparam C2 Constraint  functor 2
 * @tparam C3 Constraint  functor 3
 */
template<typename C1, typename C2, typename C3>
MetaConstraint<MetaConstraint<C1, C2>, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  MetaConstraint<C1, C2> mc;
  mc.c1 = c1;
  mc.c2 = c2;
  MetaConstraint<MetaConstraint<C1, C2>, C3> mc2;
  mc2.c1 = mc;
  mc2.c2 = c3;
  return mc2;
}


int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<Node> nodes;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      graph.add_edge(nodes[t[0]], nodes[t[1]]);
      graph.add_edge(nodes[t[0]], nodes[t[2]]);
//#if 0
      // Diagonal edges: include as of HW2 #2
      graph.add_edge(nodes[t[0]], nodes[t[3]]);
      graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
      graph.add_edge(nodes[t[1]], nodes[t[3]]);
      graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  
  // initialize Const c
  c = 1/double(graph.num_nodes());
  
  // initialize velocity to 0, and initialize mass
  NodeData nd = NodeData();
  nd.velocity = Point(0,0,0);
  nd.mass = 1/double(graph.num_nodes());
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value() = nd;
  }
  
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
   double len = (*it).length();
   (*it).value() = len;
  }
  
  
  // Construct Forces/Constraints

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0.0;
  double t_end =5; //= 5.0;
  
  auto mc = make_combined_constraint(SpecificConstraint(), PlaneConstraint(), RemoveSphereConstraint());
  
  auto mf = make_combined_force(GravityForce(), MassSpringForce(), DampingForce());
  for (double t = t_start; t < t_end; t += dt) {
    symp_euler_step(graph, t, dt, mf);
    
    mc(graph, t);
    
    // Clear the viewer 's nodes and edges
    viewer.clear();
    node_map.clear();
    
    // Update viewer with nodes ' new positions and new edges
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map );
    
    viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map );

    // Update viewer with nodes' new positions
    //viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
