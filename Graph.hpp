#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <math.h>

#include "CS207/Util.hpp"
#include "Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E = char>
class Graph {
 private:
  
  struct internal_node;  //for storing nodes and their neighbors in graph

 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Type of this graph. */
  typedef Graph graph_type;
  
  /** Synonym for the type of value stored in Node */
  typedef V node_value_type;

  /** Synonym for the type of value stored in Edge */
  typedef E edge_value_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Construct an empty graph. */
  Graph() {
    size_ = 0;
    edgesize_ = 0;
    
  }
  /** Default destructor */
  ~Graph() = default;

  /////////////
  // General //
  /////////////

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();  //clear internal nodes
    edges.clear();  //clear internal edges
    i2u.clear();    //clear all nodes
    //set sizes to be 0
    size_ = 0;
    edgesize_ = 0;
  }

  /////////////////
  // GRAPH NODES //
  /////////////////

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node:private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().point;
    }
    /** Return a position that can be modified */
    Point& position(){
      return fetch().point;
    }

    /** Return this node's index, a number in the range [0, graph_size()). */
    size_type index() const {
      return graph_->nodes[uid_].idx;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& x) const {
      if (x.graph_ == graph_ && x.uid_ == uid_)
        return true;
      //(void) x;          // Quiet compiler warning
      return false;
    }

    /** Test whether this node is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& x) const {
      if (graph_ == x.graph_)
        return (uid_ < x.uid_);
      else
        return (graph_ < x.graph_);
    }

    /** Return the value associated with this Node
     *  For rvalue operations
     */
    node_value_type& value(){
      return fetch().value;
    }

    /** Return the value associated with this Node
     *  For lvalue operations
     */
    const node_value_type& value() const {
      return fetch().value;
    }

    /** Return the number of incident edges of this node
     *  Incident are edges spawned by this node
     */
    size_type degree() const{
      return fetch().neighbors.size();
    }
    /* Return an iterator that points to the first incident edge of this node
     * If degree()==0, the returned iterator value shall not be dereferenced
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, uid_, 0);
    }
    /*  Returns an iterator referring to the past-the-end incident edge
       *  of this node
     *The past-the-end incident edge is the theoretical incident edge that 
       *  would follow the last edge. It does not point to any element, and 
       *  thus shall not be dereferenced.
     *If degree()==0, this function returns the same as edge_begin().
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, uid_, degree());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    graph_type* graph_;  //pointer to the associated graph
    size_type uid_;      //uid of the node
    // a private constructor for graph to construct a Node instance
    Node(const graph_type* graph, size_type uid)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid){
      }
    /**
     * Fetch internal node stored in graph
     * @post return.uid_ = uid_
     */
    internal_node& fetch() const {
      if (uid_ >= graph_->nodes.size())
        assert(false);
      return graph_->nodes[uid_];
    }
  };

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value_in The value associated with this node
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value_in = node_value_type()) {
    internal_node new_node;
    new_node.uid = nodes.size();
    new_node.point = position;
    new_node.value = value_in;
    new_node.idx = i2u.size();
    nodes.push_back(new_node);
    i2u.push_back(new_node.uid);
    size_ ++;
    return Node(this, nodes.size()-1); 
  }


  /** Remove a node from the graph.
   * @param[in] n Node to be removed
   * @return 1 if old has_node(n), 0 otherwise
   *
   * @post new size() == old size() - result.
   *
   * Can invalidate outstanding iterators. 
   * If old has_node(@a n), then @a n becomes invalid, as do any
   * other Node objects equal to @a n. All other Node objects remain valid.
   *
   * Complexity: Linear in size()+the degree of its neighbors.
   */
size_type remove_node(const Node& n){
  if (!has_node(n))
    return 0;
  //clean up adjacency list
  for (size_type i : nodes[n.uid_].neighbors){
    size_type pos = 0;
    //locate the position of n.uid_ in the adjacency list
    while(nodes[i].neighbors[pos] != n.uid_)
      ++pos;
    //delete uid_ in adj list
    nodes[i].neighbors.erase(nodes[i].neighbors.begin()+pos);
    nodes[i].edgevalues.erase(nodes[i].edgevalues.begin()+pos);
  }
  //update edge numbers
  edgesize_ -= nodes[n.uid_].neighbors.size();
  nodes[n.uid_].neighbors.clear();
  nodes[n.uid_].edgevalues.clear();
  //update index for valid nodes
  for (size_type i = nodes[n.uid_].idx+1; i < i2u.size(); ++i)
  {
    --nodes[i2u[i]].idx;
  }
  //update i2u
  i2u.erase(i2u.begin()+nodes[n.uid_].idx);
  --size_;
  return 1;
}
/** Remove a node from the graph.
   * @param[in] n_it Valid node_iterator
   * @return node_iterator that points to its next element if remove succeeds
       *itself if remove fails
   *
   * @post if remove succeeds new size() == old size() - 1.
   *
   * If old has_node(@a *n_it), then @a *n_it becomes invalid, as do any
   * other Node objects equal to @a *n_it. All other Node objects remain valid.
   *
   * Complexity: Linear in size()+the degree of its neighbors.
   */
node_iterator remove_node(node_iterator n_it){
  remove_node(*n_it);
  return n_it;
}

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.uid_ >= nodes.size())
      return false;
    if (i2u[nodes[n.uid_].idx] == n.uid_)
    {
      return true;
    }
    return false;
    
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i2u[i]);
  }


  template<typename G>
  void setFlag(G g){
    for (auto it = node_begin(); it != node_end(); ++it){
      if (g(*it) >-2)
        (*it).value().onBoundary = true;
      else
        (*it).value().onBoundary = false;
    }
  }

  /////////////////
  // GRAPH EDGES //
  /////////////////

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // fetch the node index stored in the graph first
      // then fetch the node from graph
      return graph_->node(graph_->nodes[node1_].idx);      
       
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // fetch the node index stored in the graph first
      // then fetch the node from graph
      return graph_->node(graph_->nodes[node2_].idx);
    }

    double length() const{
      Point diff = node1().position() - node2().position();
      return norm(diff);
    }



    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& x) const {
      if (graph_ == x.graph_ && node1_ == x.node1_ && node2_ == x.node2_)
        return true;
      else if (graph_ == x.graph_ && node1_ == x.node2_ && node2_ == x.node1_)
        return true;
      return false;
    }

    /** Test whether this edge is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The edge ordering relation must obey trichotomy: For any two edges x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Edge& x) const {
      if (graph_ == x.graph_)
        return (node1_ < x.node1_);
      else
        return (graph_ < x.graph_);
    }
    
    //return the value associated with edge for lvalue operations
    edge_value_type& value(){
      return graph_->edges[uid_].value;
    }

    //return the value associated with edge for rvalue operations
    const edge_value_type& value() const{
      return graph_->edges[uid_].value;
    }

    size_type index() const{
      return uid_;
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;  //pointer to the associated graph
    size_type node1_;      //the uid of node 1  
    size_type node2_;      //the uid of node 2
    size_type x_;  //the position in adjacency list 
    size_type y_;  //the position in adjacency list
    size_type uid_; //the unique id of the edge
    
    //private constructor for graph to construct edge instance
    Edge(const graph_type* graph, size_type node1, size_type node2, size_type uid)
      : graph_(const_cast<graph_type*>(graph)), node1_(node1), node2_(node2), uid_(uid){
      }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edgesize_;
  }
  /** Add an edge to the graph
   *@param[in] a, b are both valid nodes in the graph
   *@return The edge that has been added
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    size_type node1_uid = a.uid_;
    size_type node2_uid = b.uid_;
    //check if edge exists
    if (has_edge(a, b)){
      size_type i = 0;
      while(nodes[node1_uid].neighbors[i] != node2_uid)
        ++i;
      Edge result = Edge(this, node1_uid, node2_uid, nodes[node1_uid].edgevalues[i]);
      result.x_ = node1_uid;
      result.y_ = i;
      return result;
    }
    //if not, add a new edge
    nodes[node1_uid].neighbors.push_back(node2_uid);
    nodes[node1_uid].edgevalues.push_back(edges.size());
    nodes[node2_uid].neighbors.push_back(node1_uid);
    nodes[node2_uid].edgevalues.push_back(edges.size());
    internal_edge tempedge;
    tempedge.node1 = node1_uid;
    tempedge.node2 = node2_uid;
    edges.push_back(tempedge);
    //update edge size
    edgesize_++;
    Edge result = Edge(this, node1_uid, node2_uid, edges.size()-1);
    result.x_ = node1_uid;
    result.y_ = nodes[node1_uid].neighbors.size()-1;
    return result;
  }
  /** Remove an edge from the graph.
   * @param[in] n1, n2 Nodes of the edge to be removed
   * @return 1 if old has_edge(n1, n2), 0 otherwise
   *
   * @post new num_edges() == old num_edges() - result.
   *
   * Can invalidate outstanding iterators. 
   * If old has_edge(@a n1, @a n2), then @a e becomes invalid, as do any
   * other Edge objects equal to @a e. All other Edge objects remain valid.
   *
   * Complexity: less than degree().
   */
  size_type remove_edge(const Node& n1, const Node& n2)
  {
    if (!has_edge(n1, n2))
      return 0;
    size_type uid1 = n1.uid_;
    size_type uid2 = n2.uid_;
    size_type pos = 0;
    while(nodes[uid1].neighbors[pos] != uid2)
      ++pos;
    nodes[uid1].neighbors.erase(nodes[uid1].neighbors.begin()+pos);
    nodes[uid1].edgevalues.erase(nodes[uid1].edgevalues.begin()+pos);
    pos = 0;
    while(nodes[uid2].neighbors[pos] != uid1)
      ++pos;
    nodes[uid2].neighbors.erase(nodes[uid2].neighbors.begin()+pos);
    nodes[uid2].edgevalues.erase(nodes[uid2].edgevalues.begin()+pos);
    --edgesize_;
    return 1;
  } 
  /** Remove an edge from the graph.
   * @param[in] e Edge to be removed
   * @return 1 if old has_edge(e.node1(), e.node2()), 0 otherwise
   *
   * @post new num_edges() == old num_edges() - result.
   *
   * Can invalidate outstanding iterators. 
   * If old has_edge(@a e.node1(), @a e.node2()), then @a e becomes invalid, as do any
   * other Edge objects equal to @a e. All other Edge objects remain valid.
   *
   * Complexity: less than degree().
   */
  size_type remove_edge(const Edge& e)
  {
    return remove_edge(e.node1(), e.node1());
  }


  /** Remove an edge from the graph.
   * @param[in] e_it edge_iterator that points to the edge to be removed
   * @return edge_iterator that points to the next edge if remove succeeds;
       * edge_iterator that points to itself if otherwise
   *
   * @post new num_edges() == old num_edges() - result if remove succeeds
   *
   * Can invalidate outstanding iterators. 
   * If old has_edge(@a (*e_it).node1(), @a (*e_it).node2()), then @a e becomes invalid, as do any
   * other Edge objects equal to @a *e_it. All other Edge objects remain valid.
   *
   * Complexity: less than degree().
   */
  edge_iterator remove_edge(edge_iterator e_it)
  {
    remove_edge(*e_it);
    //if the last element in a row of adjacency list is removed
    //call ++ to go to the next valid edge
    if (e_it.nbidx_ < nodes[e_it.node_].neighbors.size())
    {
      return e_it;
    }
    else
      return ++e_it;
  }




  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type node1_uid = a.uid_;
    size_type node2_uid = b.uid_;
    //find Node b in the neighbors of Node a
    auto it = std::find(nodes[node1_uid].neighbors.begin(), nodes[node1_uid].neighbors.end(), node2_uid);
    if (it == nodes[node1_uid].neighbors.end())
      return false;
    return true;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  
  Edge edge(size_type i) const {
    return Edge(this, edges[i].node1, edges[i].node2, i);
  }


  

  ///////////////
  // Iterators //
  ///////////////

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private totally_ordered<NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    
    //return the Node pointed by the interator
    Node operator*() const{
      return graph_->node(idx_);
    }

    /*return the NodeIterator that points to the
     *next node
     *@pre uid_ < graph_->size()
     *@post (*result).index() = uid_ + 1
     */
    NodeIterator& operator++(){
      ++idx_;
      return *this;
    }
    /*Test whether this node_iterator is the same as @a x
     * two node_iterators are the same if they point to the 
       * same Node object 
     */
    bool operator==(const NodeIterator& x) const{
      if (graph_ == x.graph_ && idx_ == x.idx_)
        return true;
      return false;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type idx_;
    
 
    //private constructor for graph to construct NodeIterator instance
    NodeIterator(const graph_type* graph, size_type idx)
      : graph_(const_cast<graph_type*>(graph)), idx_(idx){
      }
  };

  // Supply definitions AND SPECIFICATIONS for:
  /* Return the node_iterator that points to the first node of the graph
   * If size()==0, the returned iterator value shall not be dereferenced
   */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  

  /*  Returns an iterator referring to the past-the-end node
   *The past-the-end node is the theoretical node that 
     *  would follow the last node. It does not point to any element, and 
     *  thus shall not be dereferenced.
   *If size()==0, this function returns the same as node_begin().
   */
  node_iterator node_end() const{
    return NodeIterator(this, i2u.size());
  }

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private totally_ordered<EdgeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // Supply definitions AND SPECIFICATIONS for:

    //return the Edge pointed by the interator
    Edge operator*() const{
      Edge edge = Edge(graph_, node_, graph_->nodes[node_].neighbors[nbidx_]
        , graph_->nodes[node_].edgevalues[nbidx_]);
      edge.x_ = node_;
      edge.y_ = nbidx_;
      return edge;
    }

    /*return the EdgeIterator that points to the next edge
     *@pre has_edge((*self).node1().position(), (*self),node2().position) == true
     *@post if (result !=edge_end()) (*result).node1().index() < (*result).node2().index()
     */
    //RI g_ != null
    //   uid1_ =< g->adj_[uid1_].size()
    //   it_<g->adj[uid1_].end()
    // END: uid1_ == g->adj_.size()
    //   uid1_ < *it_ 
    EdgeIterator& operator++(){
      //cycles until next edge is found or reaches the end
      assert(node_ < graph_->nodes.size()); 
      do{
        if (nbidx_ < (graph_->nodes[node_].neighbors.size()-1)){
          nbidx_++;
        }
        else{
          nbidx_ = 0;
          node_++;
          while(node_ != graph_->nodes.size() && graph_->nodes[node_].neighbors.size() == 0)
           node_++;
        }
      } while (node_ != graph_->nodes.size() && graph_->nodes[node_].neighbors[nbidx_] < node_); 
      return *this;
    }
    /*Test whether this edge_iterator is the same as @a x
     * two edge_iterators are the same if they point to the 
       * same Edge object 
     */
    bool operator==(const EdgeIterator& x) const{
      return (graph_ == x.graph_ && node_ == x.node_ && nbidx_ == x.nbidx_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type node_;
    size_type nbidx_;
    //private constructor for graph to construct NodeIterator instance
    EdgeIterator(const graph_type* graph, size_type node, size_type nbidx)
      : graph_(const_cast<graph_type*>(graph)), node_(node), nbidx_(nbidx){
      }
    
  };

  // Supply definitions AND SPECIFICATIONS for:
  /* Return the edge_iterator that points to the first edge of the graph
   * If num_edges()==0, the returned iterator value shall not be dereferenced
   */
  edge_iterator edge_begin() const{
    unsigned int i = 0;
    while (i < nodes.size()){
      if (nodes[i].neighbors.size() != 0)
        return EdgeIterator(this, i, 0);
    }
    return EdgeIterator(this, size(), 0);
  }

  /*  Returns an iterator referring to the past-the-end edge
   *The past-the-end edge is the theoretical edge that 
     *  does not point to any element, and 
     *  thus shall not be dereferenced.
   *If size()==0, this function returns the same as node_begin()
   *@post result.node_ = size()
   *@post result.idx_ = 0
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this, nodes.size(), 0);
  }


  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // Supply definitions AND SPECIFICATIONS for:
    //return the Edge pointed by the incident_iterator
    //@post result.node1().index() = uid_;  result.node2().index()=idx_
    Edge operator*() const{
      Edge edge = Edge(graph_, uid_, graph_->nodes[uid_].neighbors[idx_],
        graph_->nodes[uid_].edgevalues[idx_]);
      edge.x_ = uid_;
      edge.y_ = idx_;
      return edge;
    }
    /*return the Incident Iterator that points to the
     *next incident edge
     *@pre idx_<graph_->nodes[uid_].neighbors.size()
     *@post result.idx_ = idx_ + 1
     *@post if (result != edge_end()) (*result).node1().index() == idx_
     */
    IncidentIterator& operator++(){  
      idx_++;
      return *this;
    }
    /*Test whether this incident_iterator is the same as @a x
     * two incident_iterators are the same if they point to the 
       * same Edge object, and spawned by the same Node 
     */
    bool operator==(const IncidentIterator& x) const{
      return (graph_ == x.graph_ && uid_ == x.uid_ && idx_ == x.idx_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type uid_;    //the uid of the Node that spawns the incident_iterator
    size_type idx_;    //the index of the edge pointing to
    //private constructor for graph to construct NodeIterator instance
    IncidentIterator(const graph_type* graph, size_type uid, size_type idx)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid), idx_(idx){
      }
  };

 private:
  
  struct internal_node {
    Point point;
    size_type uid;  //uid of a node
    size_type idx;  //index of a node
    std::vector<size_type> neighbors;   //the uid of adjacent nodes
    node_value_type value;        //the additional value stored in nodes
    //the index that points to an element in edges
    //this vector serves as storing the unique id of an edge
    std::vector<size_type> edgevalues;  
  };

  struct internal_edge{
    edge_value_type value;
    size_type node1;
    size_type node2;
  };

  
  std::vector<internal_node> nodes;  //store node information in graph
  std::vector<size_type> i2u;   //convert the uid to index for nodes
  size_type size_;
  size_type edgesize_;
  //store edge information in graph, the position is the unique id of an edge
  std::vector<internal_edge> edges; 

};

#endif
