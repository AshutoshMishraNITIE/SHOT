//=============================================================================
/*! _zcovector*zrovector operator */
inline _zgematrix operator*(const _zcovector& covec, const zrovector& rovec)
{VERBOSE_REPORT;
  zgematrix newmat(covec.l, rovec.l);
  for(long i=0; i<newmat.m; i++){
    for(long j=0; j<newmat.n; j++){
      newmat(i,j) =covec(i)*rovec(j);
    }
  }
  
  covec.destroy();
  return _(newmat);
}
