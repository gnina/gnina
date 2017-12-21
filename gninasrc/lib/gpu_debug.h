#pragma once
#include <stdio.h>

#ifdef __CUDA_ARCH__
#define pretty_print_array(_arr, _num_elem, arr_name, elem_format)  \
    ({                                                              \
        auto arr = _arr;                                            \
        auto num_elem = _num_elem;                                  \
        if(threadIdx.x == 0){                                       \
            printf(arr_name ":[");                                  \
            for(size_t i = 0; i < num_elem; i++){                   \
                printf(elem_format ", ", arr[i]);                   \
            }                                                       \
            printf("]\n");                                          \
        }                                                           \
    })
#else
#define pretty_print_array(_arr, _num_elem, arr_name, elem_format)  \
    ({                                                              \
        auto arr = _arr;                                            \
        auto num_elem = _num_elem;                                  \
        printf(arr_name ":[");                                      \
        for(size_t i = 0; i < num_elem; i++){                       \
            printf(elem_format ", ", arr[i]);                       \
        }                                                           \
        printf("]\n");                                              \
    })
#endif

#ifdef __CUDA_ARCH__
#define pretty_print_vec_array(_arr, _num_vecs, arr_name)   \
    ({                                                      \
        vec* arr = _arr;                                    \
        sz num_vecs = _num_vecs;                            \
        if(threadIdx.x == 0){                               \
            printf(arr_name ":[");                          \
            for(size_t i = 0; i < num_vecs; i++){           \
                printf("    %lu-[%f %f %f]\n  ",             \
                       i,                                   \
                       arr[i][0],                           \
                       arr[i][1],                           \
                       arr[i][2]);                          \
            }                                               \
            printf("]\n");                                  \
        }                                                   \
    })
#else
#define pretty_print_vec_array(_arr, _num_vecs, arr_name)   \
    ({                                                      \
        vec* arr = _arr;                                    \
        sz num_vecs = _num_vecs;                            \
        printf(arr_name ":[");                              \
        for(size_t i = 0; i < num_vecs; i++){               \
            printf("    %lu-[%f %f %f]\n  ",                 \
                   i,                                       \
                   arr[i][0],                               \
                   arr[i][1],                               \
                   arr[i][2]);                              \
        }                                                   \
        printf("]\n");                                      \
    })
#endif
