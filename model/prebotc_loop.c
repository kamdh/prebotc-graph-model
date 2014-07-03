//// main routine
int num_vertices = Nvertex_types[0];
int num_edges = Nedge_list[0];
int offset = num_vertices*num_eqns_per_vertex;
int i,j,k;
//// vertex variables
for (i=0; i<num_vertices; i++) {
  double EL_cur, gP_cur, gCaN_cur, I_L, I_Na, I_K,
    I_NaP, I_CaN, I_pump, I_syn;
  int type_idx;
  j = i*num_eqns_per_vertex; // index of first (V) variable
  type_idx = vertex_types(i); // neuron type
  //// type-dependent values
  EL_cur = EL_array(type_idx); // leak potential
  gP_cur = gP_array(type_idx); // Na p conductance
  gCaN_cur = gCaN_array(type_idx); // Ca n conductance
  //// calculate currents
  // I_L(V)
  I_L = gL * ( y(j) - EL_cur );
  // I_Na(V, h, m)
  I_Na = gNa * pow(y(j+1),3) * y(j+2) * (y(j) - ENa);
  // I_K(V, n)
  I_K = gK * pow(y(j+3),4) *(y(j) - EK);
  // I_NaP(V, hp)
  I_NaP = gP_cur * infHN( Qmp, simp, y(j) ) * y(j+4) * (y(j) - ENa);
  // I_CaN(V, Ca)
  I_CaN = gCaN_cur * (y(j) - ECaN) * infHN( kCAN, siCAN, y(j+5) );
  // I_pump(Na)
  I_pump = 0.2 * ( phi( y(j+6), kNa ) - phi( Nab, kNa ) );
  //// calculate synaptic current I_syn
  double gSynTot = 0.0;
  double synVarMean = 0.0;
  for (k=0; k < (int)in_degrees(i); k++) {
    gSynTot += y(offset + (int) in_edges(i, k)) *
      edge_list( (int) in_edges(i, k), 2 );
    synVarMean += y(offset + (int) in_edges(i, k));
  }
  synVarMean = synVarMean / in_degrees(i);
  // Ohm's law
  I_syn = gSynTot * (y(j) - Esyn);
  //// set the derivatives
  // voltage V
  dydt(j) = -( Iapp + 
               I_L +
               I_Na + 
               I_K + 
               I_NaP + 
               I_CaN + 
               I_pump +
               I_syn ) / Cm;
  // Na m
  dydt(j+1) = ( infHN(Qm,sim,y(j)) - y(j+1) ) / tau(taum, Qm, sim, y(j));
  // Na h
  dydt(j+2) = ( infHN(Qh,sih,y(j)) - y(j+2) ) / tau(tauh, Qh, sih,y(j));
  // K n
  dydt(j+3) = ( infHN(Qn,siN,y(j)) - y(j+3) ) / tau(taun, Qn, siN,y(j));
  // hp Nap
  dydt(j+4) = ehp * ( infHN( Qhp, sihp, y(j) ) - y(j+4) ) /
    tau( tauhp, Qhp, sihp, y(j) );
  // Ca Can
  dydt(j+5) = eCa * ( kIP3 * synVarMean -
                      kCa * ( y(j+5) - Cab ) );
  // Na pump
  dydt(j+6) = alpha * ( - I_CaN - I_pump );
 }
// edge variables
for (i=0; i<num_edges; i++) {
  double v_source, v_target;
  v_source = y( (int)edge_list(i, 0) * num_eqns_per_vertex );
  //v_target = edge_list(i, 1);
  j = offset + i;
  dydt(j) = ( (1 - y(j)) *
              infHN( Qs, sis, v_source )  -
              ks * y(j) ) / taus;
}
