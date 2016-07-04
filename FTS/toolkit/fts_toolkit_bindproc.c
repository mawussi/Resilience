#ifdef HAVE_LIBHWLOC
#include <hwloc.h>

/*> bind an MPI process to a core.
**  Each MPI processes are separated by <topo_t> cores.
** 
** param [in] my_mpirank the mpirank of the process
** param [in] topo_n     the number of nodes
** param [in] topo_c     the number of cores per node
** param [in] topo_t     the number of threads per MPI process
** param [in] topo_p     the number of MPI processes
 */
int fts_toolkit_bind_proc_
(int * my_mpirank, int *topo_n, int * topo_c, int * topo_t, int * topo_p)
{
  hwloc_topology_t topology; /* Topology object */
  hwloc_obj_t      obj;      /* Hwloc object    */
  hwloc_cpuset_t   cpuset;   /* HwLoc cpuset    */
  int cpu;                   /* cpu rank to be computed */
  
    /* Allocate and initialize topology object.  */
    hwloc_topology_init(&topology);

    /* Perform the topology detection.  */
    hwloc_topology_load(topology);

    printf("topo %i %i %i %i \n",  *topo_n,  *topo_c,  *topo_t,  *topo_p);
	
    cpu = (*my_mpirank % (*topo_c / *topo_t)) * *topo_t;
    printf("%i-%i: Hello world\n", *my_mpirank, cpu);

    /* Get last one.  */
    obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, cpu);
    if (!obj)
      return 0;

    /* Get a copy of its cpuset that we may modify.  */
    cpuset = hwloc_cpuset_dup(obj->cpuset);

    /* Get only one logical processor (in case the core is SMT/hyperthreaded).  */
    hwloc_cpuset_singlify(cpuset);

    /* And try to bind ourself there.  */
    if (hwloc_set_cpubind(topology, cpuset, HWLOC_CPUBIND_THREAD)) {
      char *str = NULL;
      hwloc_cpuset_asprintf(&str, obj->cpuset);
      printf("Couldn't bind to cpuset %s\n", str);
      free(str);
    }

    /* Get the number at Proc level */
    cpu = obj->children[0]->os_index;

    /* Free our cpuset copy */
    hwloc_cpuset_free(cpuset);

    /* Destroy topology object.  */
    hwloc_topology_destroy(topology);

    /* return success */
    return 0;
}

#else
/* dummy function */
int fts_toolkit_bind_proc_
(int * my_mpirank, int *topo_n, int * topo_c, int * topo_t, int * topo_p)
{ return 0; }
#endif
