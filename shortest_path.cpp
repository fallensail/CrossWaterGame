/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>
#include <queue>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator
    : private equality_comparable<filter_iterator<Pred,It>> {
 public:
  // Get all of the iterator traits and make them our own
  typedef typename std::iterator_traits<It>::value_type        value_type;
  typedef typename std::iterator_traits<It>::pointer           pointer;
  typedef typename std::iterator_traits<It>::reference         reference;
  typedef typename std::iterator_traits<It>::difference_type   difference_type;
  typedef typename std::input_iterator_tag                     iterator_category;

  typedef filter_iterator<Pred,It> self_type;

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {
    // HW1 #4: YOUR CODE HERE
  }

  // HW1 #4: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return a new edge object filter_iterator itself pointing to
   *
   * @pre (*this) valid
   *
   * Complexity: O(1).
   */
  value_type operator*() const {
    return *it_;
  }
  /** Return a new edge object filter_iterator itself pointing to
   */
  self_type& operator++() {
    do {
      ++it_;
    } while (it_ != end_ && !p_(*it_));
    return *this;
  }
  /** Test whether this IncidentIterator and @a st are equal.
   *
   * Equal IncidentIterator have the same graph and the same iterator position and end position.
   */
  bool operator==(const self_type& st) const {
    return ((st.it_ == it_) && (st.end_ == end_));
  }

 private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

// HW1 #4: YOUR CODE HERE
// Specify and write an interesting predicate on the nodes.
// Explain what your predicate is intended to do and test it.
// If you'd like you may create new nodes and tets files.

// this predicate will only get the largest cross section
struct MyNodePredicate {
  template <typename Node>
  bool operator() (const Node& n) const {
    return (n.position().x) < 0.01 && (n.position().x > -0.01); //dist < 0.000001;
  }
};

/** Test predicate for HW1 #4 */
struct SlicePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
  //std::cout<<"SlicePredicate called: " << (n.position().x < 0) << std::endl;
    return n.position().x < 0;
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
    //hw1 5
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
  // HW1 #4: YOUR CODE HERE
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

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();
  
  // HW1 #4: YOUR CODE HERE
  
  // Use shortest_path_lengths to set the node values to the path lengths
  int max = shortest_path_lengths(graph, Point(-1,0,1));
  auto node_map = viewer.empty_node_map(graph);
  
  // Construct a Color functor and view with the SDLViewer
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
  
  ColorFn cf = ColorFn(max);
  
  viewer.add_nodes(graph.node_begin(), graph.node_end(), cf, node_map);
  //MyNodePredicate predicate;
  //viewer.add_nodes(make_filtered(graph.node_begin(), graph.node_end(), predicate), make_filtered(graph.node_end(), graph.node_end(), predicate), cf, node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();
  
  return 0;
}
