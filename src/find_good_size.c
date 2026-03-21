/* < find a number nn such that n<nn, nn=MIN(2^m, 3/4*2^m, 5/8*2^m >
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/

int find_good_size(int n)
{
  int m, p, nn;

  m = 1;
  while(m<=n) m *= 2;//m>=n

  nn = m;  
  /* p = 7*m/8; */
  /* if(p>n && p<nn) nn = p; */
  p = 3*m/4;
  if(p>n && p<nn) nn = p;
  p = 5*m/8;
  if(p>n && p<nn) nn = p;
    
  return nn;
}
