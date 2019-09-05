#include <iomanip>
#include "nto_analysis.h"
#include <libwfa/soc.h>

namespace libwfa {

using namespace arma;

nto_analysis::nto_analysis(const mat &s, const ab_matrix &c,
                           const h_so& h_so1e, const h_so& h_somf, 
                           const arma::mat &tdm) {
    initialize_proper(s, c, tdm);
    //std::cout << "SOC will be soon.." << std::endl;
    std::vector<soc> vec_so1e;
    std::vector<soc> vec_somf;
    const double cm1 = 219474.63;

    //  Reciprocal speed of light in a.u.
    const double rc = 2.1876912633e6/299792458.0;
    //  Breit-Pauli equation prefactor: 1/(2c^2)
    const double prefac = 0.5 * rc * rc;

    arma::Col<std::complex<double> > so1e_int_lm1_sp1;
    arma::Col<std::complex<double> > so1e_int_l0_s0;
    arma::Col<std::complex<double> > so1e_int_lp1_sm1;


    {
    arma::Mat<std::complex<double> > tmp_mat1 = m_nto[0]->get_coeff().t()*h_so1e.lm1_sp1*m_nto[1]->get_coeff();
    so1e_int_lm1_sp1 = tmp_mat1.diag();
    arma::Mat<std::complex<double> > tmp_mat2 = m_nto[0]->get_coeff().t()*h_so1e.l0_s0*m_nto[1]->get_coeff();
    so1e_int_l0_s0 = tmp_mat2.diag();
    arma::Mat<std::complex<double> > tmp_mat3 = m_nto[0]->get_coeff().t()*h_so1e.lp1_sm1*m_nto[1]->get_coeff();
    so1e_int_lp1_sm1 = tmp_mat3.diag();
    }
    vec_so1e.resize(so1e_int_lm1_sp1.size());
    for(size_t i = 0; i < vec_so1e.size(); i++) {
        vec_so1e[i].lm1_sp1  =  so1e_int_lm1_sp1(i);
        vec_so1e[i].l0_s0    =  so1e_int_l0_s0(i);
        vec_so1e[i].lp1_sm1  =  so1e_int_lp1_sm1(i);
    }
    m_nto[0]->get_so1e() = vec_so1e;

    arma::Col<std::complex<double> > somf_int_lm1_sp1;
    arma::Col<std::complex<double> > somf_int_l0_s0;
    arma::Col<std::complex<double> > somf_int_lp1_sm1;

    {
    arma::Mat<std::complex<double> > tmp_mat1 = m_nto[0]->get_coeff().t()*h_somf.lm1_sp1*m_nto[1]->get_coeff();
    somf_int_lm1_sp1 = tmp_mat1.diag();
    arma::Mat<std::complex<double> > tmp_mat2 = m_nto[0]->get_coeff().t()*h_somf.l0_s0*m_nto[1]->get_coeff();
    somf_int_l0_s0 = tmp_mat2.diag();
    arma::Mat<std::complex<double> > tmp_mat3 = m_nto[0]->get_coeff().t()*h_somf.lp1_sm1*m_nto[1]->get_coeff();
    somf_int_lp1_sm1 = tmp_mat3.diag();
    }
    vec_somf.resize(somf_int_lm1_sp1.size());
    for(size_t i = 0; i < vec_somf.size(); i++) {
        vec_somf[i].lm1_sp1  =  somf_int_lm1_sp1(i);
        vec_somf[i].l0_s0    =  somf_int_l0_s0(i);
        vec_somf[i].lp1_sm1  =  somf_int_lp1_sm1(i);
    }
    m_nto[0]->get_somf() = vec_somf;
}

nto_analysis::nto_analysis(const mat &s, const ab_matrix &c,
    const ab_matrix &tdm) {
  
#if 0
  //Leave the old code just in case someone wants to revert
  ab_matrix edm, hdm;
  form_eh(s, tdm, edm, hdm);
  initialize(s, c, edm, hdm);
#else
    //TDM is hole-particle, i.e., \gamma^{IF}
    initialize_proper(s, c, tdm);
#endif
}

nto_analysis::nto_analysis(const mat &s, const ab_matrix &c, const ab_matrix &fock,
			   const ab_matrix &tdm, bool use_fock) {
  
  initialize_proper(s, c, tdm);

  if (use_fock) {  //Now intitalize energies

    arma::vec ene;
    arma::mat tmp_mat;

#if 0
    //Test that Fock matrix was read coorectly
    arma::mat cs;
    cs=c.alpha();
    tmp_mat=cs.t()*fock.alpha()*cs;
    ene=tmp_mat.diag();
    std::cout<< "Diag: c+Fc" << std::endl;
    ene.print();
    //Check if normalized. NTOs in cols...
    tmp_mat=m_nto[0]->get_coeff().t()*s*m_nto[0]->get_coeff();
    ene=tmp_mat.diag();
    std::cout<< "U+SU" << std::endl;
    ene=tmp_mat.diag();
    ene.print();
#endif

    //Transform Fock matrix into hole NTO basis. NTOs stored by cols...
    //First particles, then holes
    tmp_mat=m_nto[0]->get_coeff().t()*fock.alpha()*m_nto[0]->get_coeff();
    ene=tmp_mat.diag();
    m_nto[0]->get_ene()=ene;
    
    //Now same for holes...
    tmp_mat=m_nto[1]->get_coeff().t()*fock.alpha()*m_nto[1]->get_coeff();
    ene=tmp_mat.diag();
    m_nto[1]->get_ene()=ene;
    
    if (m_nto[2]) {
      
      //do the same for beta
      tmp_mat=m_nto[2]->get_coeff().t()*fock.alpha()*m_nto[2]->get_coeff();
      ene=tmp_mat.diag();
      m_nto[2]->get_ene()=ene;
      tmp_mat=m_nto[3]->get_coeff().t()*fock.alpha()*m_nto[3]->get_coeff();
      ene=tmp_mat.diag();
      m_nto[3]->get_ene()=ene;
    }
  }
}

  void nto_analysis::analyse(std::ostream &out, size_t nnto) const {

    if (nnto == 0)  return;

    bool got_ene=m_nto[0]->got_ene();
    bool got_soc=m_nto[0]->got_soc();

    std::vector<soc>& vec_so1e = m_nto[0]->get_so1e();
    
    if(!got_soc) {
        if (m_nto[2]) {
          out << "NTOs (alpha)" << std::endl;
          if (got_ene)
            analysis(out, m_nto[0]->get_occ(), m_nto[1]->get_ene(), m_nto[0]->get_ene(), nnto);
          else
            analysis(out, m_nto[0]->get_occ(), nnto);
            
          out << "NTOs (beta)" << std::endl;
          if (got_ene)
            analysis(out, m_nto[2]->get_occ(), m_nto[3]->get_ene(), m_nto[2]->get_ene(), nnto);
          else
            analysis(out, m_nto[2]->get_occ(), nnto);
        }
        else {
          out << "NTOs" << std::endl;
          if (got_ene)
            analysis(out, m_nto[0]->get_occ()*2.0, m_nto[1]->get_ene(), m_nto[0]->get_ene(), nnto);
          else
            analysis(out, m_nto[0]->get_occ()*2.0, nnto);
        }
    }
    else {
          out << "spinless NTOs of universal triplet transition OPDM" << std::endl;
          analysis(out, m_nto[0]->get_so1e(), m_nto[0]->get_somf(), m_nto[0]->get_occ(), nnto);
    }
}

void nto_analysis::export_orbitals(orbital_printer_i &pr, double thresh) const {

    if (m_nto[2]) {

        const vec &ee_a = m_nto[0]->get_occ(),   &eh_a = m_nto[1]->get_occ();
        const mat &ce_a = m_nto[0]->get_coeff(), &ch_a = m_nto[1]->get_coeff();
        const vec &ee_b = m_nto[2]->get_occ(),   &eh_b = m_nto[2]->get_occ();
        const mat &ce_b = m_nto[3]->get_coeff(), &ch_b = m_nto[3]->get_coeff();
        orbital_data nto_a(join_cols(eh_a * -1., flipud(ee_a)),
                join_rows(ch_a, fliplr(ce_a)));
        orbital_data nto_b(join_cols(eh_b * -1., flipud(ee_b)),
                join_rows(ch_b, fliplr(ce_b)));
        orbital_selector s_a, s_b;
        build_selector(ee_a, eh_a, thresh, s_a);
        build_selector(ee_b, eh_b, thresh, s_b);

        pr.perform(orbital_type::nto, nto_a, s_a, nto_b, s_b);
    }
    else {

        const vec &ee = m_nto[0]->get_occ(), &eh = m_nto[1]->get_occ();
        const mat &ce = m_nto[0]->get_coeff(), &ch = m_nto[1]->get_coeff();
        orbital_data nto(join_cols(eh * -1., flipud(ee)),
                join_rows(ch, fliplr(ce)));
        orbital_selector s;
        build_selector(ee, eh, thresh, s);

        pr.perform(orbital_type::nto, nto, s);
    }
}


void nto_analysis::form_eh(const mat &s, const ab_matrix &tdm,
    ab_matrix &edm, ab_matrix &hdm) {

    if (tdm.is_alpha_eq_beta()) {
        edm.set_alpha_eq_beta();
        hdm.set_alpha_eq_beta();

        form_eh(s, tdm.alpha(), edm.alpha(), hdm.alpha());
    }
    else {
        edm.set_alpha_neq_beta();
        hdm.set_alpha_neq_beta();

        form_eh(s, tdm.alpha(), edm.alpha(), hdm.alpha());
        form_eh(s, tdm.beta(), edm.beta(), hdm.beta());
    }
}


void nto_analysis::initialize(const arma::mat &s, const ab_matrix &c,
        const ab_matrix &edm, const ab_matrix &hdm) {

    if (edm.is_alpha_eq_beta()) {
        m_nto[0] = new orbital_data(s, c.alpha(), edm.alpha());
        m_nto[1] = new orbital_data(s, c.alpha(), hdm.alpha());
        m_nto[2] = m_nto[3] = 0;
    }
    else {
        m_nto[0] = new orbital_data(s, c.alpha(), edm.alpha());
        m_nto[1] = new orbital_data(s, c.alpha(), hdm.alpha());
        m_nto[2] = new orbital_data(s, c.beta(), edm.beta());
        m_nto[3] = new orbital_data(s, c.beta(), hdm.beta());
    }
}


  //Compute NTOs without doing square trick, but using proper SVD procedure
void nto_analysis::initialize_proper(const arma::mat &s, const ab_matrix &c,
				     const ab_matrix &tdm) {

  //gamma_mo=C'xSx gamma_aox SxC
  arma::mat cs=s*c.alpha();
  arma::mat tdm_mo=cs.t()*tdm.alpha()*cs;
  //arma::mat tdm_mo=c.alpha().t()*s*tdm.alpha()*s*c.alpha();
  
  arma::mat alpha, beta;
  arma::vec sigma;
  //tdm_mo=alpha x sigma x beta_t
  //Can use different methods: "dc" or "std". "dc" IS faster.
  arma::svd(alpha,sigma,beta,tdm_mo,"dc");

  //libwfa wants sigma^2 and it wants them in ascending order
  size_t nsigmas=sigma.n_elem;
  arma::vec sigma_reorder(nsigmas);
  for(size_t i=0; i<nsigmas; i++) 
    sigma_reorder[i]=sigma[nsigmas-1-i]*sigma[nsigmas-1-i]; 

#if 0
  //Loss of precision occur here...
  double nrm1=arma::norm(tdm_mo);
  nrm1=nrm1*nrm1*2.0; 
  std::cout << std::endl <<
    "Omega_g= " << nrm1 << "  ||gamma_g||=" << sqrt(nrm1) << std::endl;
  double nrm2=arma::norm(sigma);
  nrm2=nrm2*nrm2*2.0; 
  std::cout << "Omega_s= " << nrm2 << "  ||gamma_s||=" << sqrt(nrm2) << std::endl;
#endif  
    
  //This now contains particle NTOs in AO basis
  tdm_mo=c.alpha()*beta;
  beta.resize(tdm_mo.n_rows,tdm_mo.n_cols);
  //Stupid step, but need to do it for libwfa
  for(size_t i=0; i<nsigmas; i++)
    beta.col(i)=tdm_mo.col(nsigmas-1-i);

  m_nto[0] = new orbital_data(sigma_reorder, beta);

  //This now contains hole NTOs in AO basis, NTOs arranged by coloumns
  tdm_mo=c.alpha()*alpha;
  alpha.resize(tdm_mo.n_rows,tdm_mo.n_cols);
  for(size_t i=0; i<nsigmas; i++)
    alpha.col(i)=tdm_mo.col(nsigmas-1-i);

  m_nto[1] = new orbital_data(sigma_reorder,alpha);

  if (tdm.is_alpha_eq_beta()) { 
    //we are done
    m_nto[2] = m_nto[3] = 0;
  }
  else {
    
    //repeat the procedure for beta part
    cs=s*c.beta();
    tdm_mo=cs.t()*tdm.beta()*cs;
    //tdm_mo=c.beta().t()*s*tdm.beta()*s*c.beta();
    
    //tdm_mo=alpha x sigma x beta_t
    arma::svd(alpha,sigma,beta,tdm_mo,"dc");
    
    //libwfa wants sigma^2 
    //libwfa wants sigma^2 and it wants them in ascending order
    for(size_t i=0; i<nsigmas; i++) 
      sigma_reorder[i]=sigma[nsigmas-1-i]*sigma[nsigmas-1-i]; 
    
    //This now contains particle NTOs in AO basis
    tdm_mo=c.beta()*beta;
    beta.resize(tdm_mo.n_rows,tdm_mo.n_cols);
    for(size_t i=0; i<nsigmas; i++)
      beta.col(i)=tdm_mo.col(nsigmas-1-i);
    m_nto[2] = new orbital_data(sigma_reorder, beta);
    
    //This now contains hole NTOs in AO basis
    tdm_mo=c.beta()*alpha;
    alpha.resize(tdm_mo.n_rows,tdm_mo.n_cols);
    for(size_t i=0; i<nsigmas; i++)
      alpha.col(i)=tdm_mo.col(nsigmas-1-i);
    m_nto[3] = new orbital_data(sigma_reorder,alpha);
  }
  
}


  //Compute NTOs without doing square trick, but using proper SVD procedure
  //This is a version for a spinless matrix for SOC
void nto_analysis::initialize_proper(const arma::mat &s, const ab_matrix &c,
				     const arma::mat &tdm) {

  //gamma_mo=C'xSx gamma_aox SxC
  arma::mat cs=s*c.alpha();
  arma::mat tdm_mo=cs.t()*tdm*cs;
  
  arma::mat alpha, beta; // PP: this notation is super-confusing. Alpha is the left vectors, beta - the right vectors
  arma::vec sigma;
  //tdm_mo=alpha x sigma x beta_t
  //Can use different methods: "dc" or "std". "dc" IS faster.
  arma::svd(alpha,sigma,beta,tdm_mo,"dc");

  //libwfa wants sigma^2 and it wants them in ascending order
  size_t nsigmas=sigma.n_elem;
  arma::vec sigma_reorder(nsigmas);
  for(size_t i=0; i<nsigmas; i++) 
    sigma_reorder[i]=sigma[nsigmas-1-i]*sigma[nsigmas-1-i]; 

  //This now contains particle NTOs in AO basis
  tdm_mo=c.alpha()*beta;
  beta.resize(tdm_mo.n_rows,tdm_mo.n_cols);
  //Stupid step, but need to do it for libwfa
  for(size_t i=0; i<nsigmas; i++)
    beta.col(i)=tdm_mo.col(nsigmas-1-i);

  m_nto[0] = new orbital_data(sigma_reorder, beta);

  //This now contains hole NTOs in AO basis, NTOs arranged by coloumns
  tdm_mo=c.alpha()*alpha;
  alpha.resize(tdm_mo.n_rows,tdm_mo.n_cols);
  for(size_t i=0; i<nsigmas; i++)
    alpha.col(i)=tdm_mo.col(nsigmas-1-i);

  m_nto[1] = new orbital_data(sigma_reorder,alpha);
    //we are done
  m_nto[2] = m_nto[3] = 0;
}
  
void nto_analysis::analysis(std::ostream &out,
    const arma::vec &e, size_t nnto) {

    out << "  Leading SVs:" << std::endl;
    out << std::setprecision(4) << std::fixed;
    out << "  ";
    for (size_t i = 0, j = e.n_rows - 1; i < nnto; i++, j--) {
        out << std::setw(9) << e(j);
    }
    out << std::endl;

    double Om = accu(e);

    out << std::setprecision(6) << std::fixed;
    out << "  Sum of SVs (Omega):            "
        << std::setw(11) << Om << std::endl;
    out << "  Participation ratio (PR_NTO):  "
        << std::setw(11) << Om * Om / dot(e, e);
    out << std::endl;

    // Entanglement values
    // Compute and print only for reasonable Omega
    if (Om < 1.5) {
        double ln2 = 0.6931471805599453;
        double SHE  = -arma::dot(e, arma::trunc_log(e));
            SHE /= ln2;
        double rSHE = -arma::dot(e/Om, arma::trunc_log(e/Om));
            rSHE /= ln2;

        out << "  Entanglement entropy (S_HE):   "
            << std::setw(11) << SHE << std::endl;
        out << "  Nr of entangled states (Z_HE): "
            << std::setw(11) << pow(2, SHE) << std::endl;
        out << "  Renormalized S_HE/Z_HE:"
            << std::setw(10) << rSHE << " /"
            << std::setw(10) << pow(2, rSHE) << std::endl;
    }
}

void nto_analysis::analysis(std::ostream &out, 
    std::vector<soc>& vec_so1e, std::vector<soc>& vec_somf,
    const arma::vec &e, size_t nnto) {

    const double cm1 = 219474.63;

    //  Reciprocal speed of light in a.u.
    const double rc = 2.1876912633e6/299792458.0;
    //  Breit-Pauli equation prefactor: 1/(2c^2)
    const double prefac = 0.5 * rc * rc;

    out << "  Leading SVs:" << std::endl;
    out << std::setprecision(4) << std::fixed;
    out << "  ";
    for (size_t i = 0, j = e.n_rows - 1; i < nnto; i++, j--) {
        out << std::setw(9) << e(j);
    }
    out << std::endl;

    // PP: it is so hard to trace these signes...
    // Phases are complicated because of the property of conjugation of spin-tensor

    out << "  1e SOC integrals on NTOs: (with prefactors, cm-1)" << std::endl;
    out << std::setprecision(4) << std::fixed;
    out << "L-";
    for (size_t i = 0, j = vec_so1e.size() - 1; i < nnto; i++, j--) {
        out << std::setw(9) << +0.5*std::conj(vec_so1e[j].lm1_sp1)*cm1*prefac;
    }
    out << std::endl;
    out << "L0";
    for (size_t i = 0, j = vec_so1e.size() - 1; i < nnto; i++, j--) {
        out << std::setw(9) << -(std::sqrt(2.0)/2.0)*std::conj(vec_so1e[j].l0_s0)*cm1*prefac;
    }
    out << std::endl;
    out << "L+";
    for (size_t i = 0, j = vec_so1e.size() - 1; i < nnto; i++, j--) {
        out << std::setw(9) << -0.5*std::conj(vec_so1e[j].lp1_sm1)*cm1*prefac;
    }
    out << std::endl;

    out << "  Mean-field SOC integrals on NTOs: (with prefactors, cm-1)" << std::endl;
    out << std::setprecision(4) << std::fixed;
    out << "L-";
    for (size_t i = 0, j = vec_somf.size() - 1; i < nnto; i++, j--) {
        out << std::setw(9) << +0.5*std::conj(vec_somf[j].lm1_sp1)*cm1*prefac;
    }
    out << std::endl;
    out << "L0";
    for (size_t i = 0, j = vec_somf.size() - 1; i < nnto; i++, j--) {
        out << std::setw(9) << -(std::sqrt(2.0)/2.0)*std::conj(vec_somf[j].l0_s0)*cm1*prefac;
    }
    out << std::endl;
    out << "L+";
    for (size_t i = 0, j = vec_somf.size() - 1; i < nnto; i++, j--) {
        out << std::setw(9) << -0.5*std::conj(vec_somf[j].lp1_sm1)*cm1*prefac;
    }
    out << std::endl;


    double Om = accu(e);

    out << std::setprecision(6) << std::fixed;
    out << "  Sum of SVs (Omega):            "
        << std::setw(11) << Om << std::endl;
    out << "  Participation ratio (PR_NTO):  "
        << std::setw(11) << Om * Om / dot(e, e);
    out << std::endl;

    // Entanglement values
    // Compute and print only for reasonable Omega
    if (Om < 1.5) {
        double ln2 = 0.6931471805599453;
        double SHE  = -arma::dot(e, arma::trunc_log(e));
            SHE /= ln2;
        double rSHE = -arma::dot(e/Om, arma::trunc_log(e/Om));
            rSHE /= ln2;

        out << "  Entanglement entropy (S_HE):   "
            << std::setw(11) << SHE << std::endl;
        out << "  Nr of entangled states (Z_HE): "
            << std::setw(11) << pow(2, SHE) << std::endl;
        out << "  Renormalized S_HE/Z_HE:"
            << std::setw(10) << rSHE << " /"
            << std::setw(10) << pow(2, rSHE) << std::endl;
    }
}


  void nto_analysis::analysis(std::ostream &out,
			      const arma::vec &e,
			      const arma::vec &hole_e,
			      const arma::vec &e_e,
			      size_t nnto) {

    
    out << "  Leading SVs and ENEs:" << std::endl;
    out << std::setprecision(4) << std::fixed;
    out << "    ";
    for (size_t i = 0, j = e.n_rows - 1; i < nnto; i++, j--)
        out << std::setw(9) << e(j);
    out << std::endl;
    out << "E(h)";
    for (size_t i = 0, j = hole_e.n_rows - 1; i < nnto; i++, j--)
      out << std::setw(9) << hole_e(j);
    out << std::endl;
    out << "E(e)";
    for (size_t i = 0, j = e_e.n_rows - 1; i < nnto; i++, j--)
      out << std::setw(9) << e_e(j);
    out << std::endl << std::endl;

    double Om = accu(e);

    out << std::setprecision(6) << std::fixed;
    out << "  Sum of SVs (Omega):            "
        << std::setw(11) << Om << std::endl;
    out << "  Participation ratio (PR_NTO):  "
        << std::setw(11) << Om * Om / dot(e, e);
    out << std::endl;

    // Entanglement values
    // Compute and print only for reasonable Omega
    if (Om < 1.5) {
        double ln2 = 0.6931471805599453;
        double SHE  = -arma::dot(e, arma::trunc_log(e));
            SHE /= ln2;
        double rSHE = -arma::dot(e/Om, arma::trunc_log(e/Om));
            rSHE /= ln2;

        out << "  Entanglement entropy (S_HE):   "
            << std::setw(11) << SHE << std::endl;
        out << "  Nr of entangled states (Z_HE): "
            << std::setw(11) << pow(2, SHE) << std::endl;
        out << "  Renormalized S_HE/Z_HE:"
            << std::setw(10) << rSHE << " /"
            << std::setw(10) << pow(2, rSHE) << std::endl;
    }
}

void nto_analysis::build_selector(const arma::vec &e, const arma::vec &h,
    double thresh, orbital_selector &sel) {

    size_t ntot = h.size() + e.size();
    uvec ph = find(h > thresh, 1), pe = find(e > thresh, 1);

    size_t nh = (ph.size() == 1 ? ph(0) : 0);
    size_t ne = (pe.size() == 1 ? pe(0) : 0);
    if (sel.n_indexes() != ntot) sel = orbital_selector(ntot);

    if (nh > 0) sel.select(true, nh, h.size(), 1);
    if (ne > 0) sel.select(false, h.size(), h.size() + (e.size() - ne), 1);
}


} // namespace libwfa


