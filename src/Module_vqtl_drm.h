/*
* C implementation of Deviation Regression Model (DRM).
* The original R code is available at https://github.com/drewmard/DRM.
* -- Benjamin Fang 20230117
*/

#ifndef MODULE_VQTL_DRM_HEAD

#define MODULE_VQTL_DRM_HEAD
#if defined MODULE_VQTL_DRM_HEAD
    #define MODULE_VQTL_DRM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" {
#endif

MODULE_VQTL_DRM_EXTERN int Module_vqtl_drm(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
