










# 3. combine diagnostics from all times in cs.dist_var (DistributionChunked)

def combineMcsAgeDiagnostic(cs,diag,metric):
    """Compute the mean age"""
    
    # name of diagnostic
    diag_name = '%s_%s'%(diag,metric)
    N_t = len(cs.N_all_t)
    N_Q = cs.dist_pr_all.nbins
    
    if diag == 'mean':
        
        N_valid = np.full(N_Q,0)
        diag_value = np.full(N_Q,0)
        
        for i_t in range(N_t):

            if cs.diag_all_t[i_t] is not None: # None should correspond to indices with corrupted data
                
                N_valid = N_valid + cs.N_valid_all_t[i_t]
                diag_value = diag_value + (cs.N_valid_all_t[i_t])*(cs.diag_all_t[i_t])
        
        # normalize sum by total number of valid points in each bin
        diag_value = diag_value/N_valid
        
    return diag_value, N_valid










mean_norm_age, N_valid = combineMcsAgeDiagnostic(cs,'mean','norm_age')


##-- define slice for valid bins

# start after last nan value
# i_min = np.where(np.isnan(mean_norm_age))[0][-1]
i_min = 7
# define slice
s_hide = slice(None,i_min)
# hide noisy/nan beginning
mean_norm_age[s_hide] = np.nan
N_valid[s_hide] = np.nan