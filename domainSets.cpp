//
//  domainSets.cpp
//  uncertainGraph
//
//  Created by dongqingxiao on 1/18/17.
//  Copyright © 2017 dongqingxiao. All rights reserved.
//

#include "domainSets.hpp"
using namespace arma;
using namespace std;

int ITER = 0;

uvec is_member(vec vec,mat mat)
{
    uvec membership(vec.n_elem);
    
    for(int i=0; i<vec.n_elem; i++)
    {
        uvec a = find(mat==vec(i));
        membership(i) = a.is_empty()?0:1;
    }
    
    return membership;
}

umat dominantSetsClustering(fmat S,int N_MAX_CLUSTER_SIZE,int N_MIN_CLUSTER_SIZE, int N_OVERLAP_SIZE)
{
    
    //TIMING INITIALIZATIONS
    struct timeval tp;
    double t_total(0), start, end, sec, usec;
    
    //Time start
    gettimeofday( &tp, NULL );
    sec = static_cast<double>( tp.tv_sec );
    usec = static_cast<double>( tp.tv_usec )/1E6;
    start = sec + usec;
    
    
    
    int Nc = S.n_rows;
    umat clusters_matrix = zeros<umat>(Nc,Nc);
    umat overlaps_matrix = zeros<umat>(Nc,Nc);
    uvec cluster_assign = zeros<uvec>(Nc);
    uvec onetoNc = conv_to< uvec >::from(linspace(0,Nc-1,Nc));
    uvec rem_list = onetoNc;
    
    int i=0;
    uvec ind;
    while(rem_list.n_elem > N_MIN_CLUSTER_SIZE)
    {
        
        uvec dominant_set;
        uvec non_dominant_set;
        dominantSetExtraction(S(rem_list,rem_list),dominant_set,non_dominant_set,N_MAX_CLUSTER_SIZE,N_MIN_CLUSTER_SIZE);
        uvec new_cluster = rem_list(dominant_set);
        ind << i;
        clusters_matrix(ind,new_cluster) = ones<umat>(1,new_cluster.n_elem);
        cluster_assign(new_cluster) = (i+1)*ones<uvec>(new_cluster.n_elem);
        //Overlap constraint - Diverse Overlapping
        if(N_OVERLAP_SIZE!=0)
        {
            int next_overlap;
            uvec next_overlaps_available, sortedind, next_overlap_uvec;
            for(int n=0;n<N_OVERLAP_SIZE;n++)
            {
                if(n==0)
                {
                    next_overlaps_available = new_cluster(find(is_member(conv_to<vec>::from(new_cluster),conv_to<mat>::from(overlaps_matrix))==0));
                    next_overlap = next_overlaps_available(next_overlaps_available.n_elem-1);
                    overlaps_matrix(i,n) = next_overlap;
                }
                else
                {
                    next_overlap_uvec << next_overlap;
                    sortedind = conv_to< uvec >::from(sort_index(S(next_overlap_uvec,new_cluster),0));//ascending
                    sortedind = sortedind(find(is_member(conv_to<vec>::from(new_cluster(sortedind)),conv_to<mat>::from(overlaps_matrix))==0));
                    next_overlap = new_cluster(sortedind(0));
                    overlaps_matrix(i,n) = next_overlap;
                }
            }
        }
        
        rem_list = rem_list(non_dominant_set);
        i++;
    }
    
    int no_groups = i;
    clusters_matrix.resize(no_groups,Nc);
    overlaps_matrix.resize(no_groups,N_OVERLAP_SIZE);
    
    //Force assignment of remaining elements
    for(int i=0;i<rem_list.n_elem;i++)
    {
        int element = rem_list(i);
        uvec element_uvec;
        element_uvec << element;
        uvec elements_in_groups = find(cluster_assign!=0);
        uvec sortedind = conv_to<uvec>::from(sort_index(S(element_uvec,elements_in_groups),1));
        int nearest_group = cluster_assign(elements_in_groups(sortedind(0)));
        clusters_matrix(nearest_group-1,element) = 1;
        cluster_assign(element) = nearest_group;
    }
    
    //Assign overlaps to both groups
    for(int i=0;i<no_groups;i++)
    {
        for(int j=0;j<N_OVERLAP_SIZE;j++)
        {
            int element = overlaps_matrix(i,j);
            uvec element_uvec;
            element_uvec << element;
            uvec elements_other_groups = find(cluster_assign!=(i+1));
            uvec sortedind = conv_to<uvec>::from(sort_index(S(element_uvec,elements_other_groups),1));
            int overlap_group = cluster_assign(elements_other_groups(sortedind(0)));
            clusters_matrix(overlap_group-1,element) = 1;
        }
    }
    cout << no_groups << " clusters found";
    
    //Time end
    gettimeofday( &tp, NULL );
    sec = static_cast<double>( tp.tv_sec );
    usec = static_cast<double>( tp.tv_usec )/1E6;
    end = sec + usec;
    t_total += (end-start);
    
    
    //Print times
    cout << "\n\nDomSet Clustering times\n--------------------------------\n";
    cout << "Clustering time: " << t_total << "\n";
    
    return clusters_matrix;
}

int dominantSetExtraction(fmat S,uvec& dominant, uvec& non_dominant, int N_MAX_CLUSTER_SIZE, int N_MIN_CLUSTER_SIZE)
{
    
    
    int Nc = S.n_rows;
    S = S - conv_to<fmat>::from(eye(Nc,Nc));
    vec X = (1.0/Nc) * ones<vec>(Nc);
    double f_0 = as_scalar(X.t() * S * X);
    double f_1 = 0;
    double f_2 = f_0;
    double epsilon = 0.001;
    int n_iter = 0;
    while ((f_2-f_1)/f_0 > epsilon)
    {
        n_iter++;
        f_1 = f_2;
        X = X % ((S.t() * X)/f_1);
        f_2 = as_scalar(X.t() * S * X);
    }
    
    
    uvec ordered_list = conv_to< uvec >::from(sort_index(X,1));//descending
    uvec gteps = find(X>10e-13);
    int cluster_size = gteps.n_elem;
    
    //Max/Min size constraints
    if (cluster_size <= N_MIN_CLUSTER_SIZE)
    {
        dominant = ordered_list.rows(0,N_MIN_CLUSTER_SIZE-1);
        non_dominant = ordered_list.rows(N_MIN_CLUSTER_SIZE,Nc-1);
    }
    else if ((cluster_size > N_MIN_CLUSTER_SIZE)  && (cluster_size <= N_MAX_CLUSTER_SIZE))
    {
        dominant = ordered_list.rows(0,cluster_size-1);
        if(cluster_size!=Nc) non_dominant = ordered_list.rows(cluster_size,Nc-1);
    }
    else
    {
        dominant = ordered_list.rows(0,N_MAX_CLUSTER_SIZE-1);
        if(cluster_size!=Nc) non_dominant = ordered_list.rows(N_MAX_CLUSTER_SIZE,Nc-1);
    }
    
    
    int corrected_cluster_size = dominant.n_elem;
    cout << "Cluster size: " << corrected_cluster_size << "	Niter = " << n_iter << endl;
    return 0;
}
