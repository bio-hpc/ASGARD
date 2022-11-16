         if( ipimd>0 ) then
            if(cnum(i).eq.0.and.cnum(j).eq.0) then
               nrg_all(1:nbead)=nrg_all(1:nbead) + ecur*nbead_inv
            else 
               if(cnum(i).ne.0) then
                  nrg_all(cnum(i)) = nrg_all(cnum(i)) + ecur
               else
                  nrg_all(cnum(j)) = nrg_all(cnum(j)) + ecur
               end if
            end if
         end if
