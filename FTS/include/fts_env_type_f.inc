  Type fts_env_t
     sequence
     Integer :: comm            ! Global MPI communicator
     Integer :: nodecomm        ! MPI communicator on the node.
     Integer :: nodeid       ! color of the node.
  End type fts_env_t
