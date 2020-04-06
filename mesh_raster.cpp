#include <iostream>
#include <stdio.h>
#include <fstream>
#include <chrono>

#include "glm/vec3.hpp"
#include "glm/vec2.hpp"
#include "glm/mat4x4.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/matrix.hpp"

#include "cow.h"
using namespace std;

static const float inchToMm = 25.4; 
enum FitResolutionGate { kFill = 0, kOverscan };

void computeScreenCoordinates( 
    const float &filmApertureWidth, 
    const float &filmApertureHeight, 
    const uint32_t &imageWidth, 
    const uint32_t &imageHeight, 
    const FitResolutionGate &fitFilm, 
    const float &nearClippingPLane, 
    const float &focalLength, 
    float &top, float &bottom, float &left, float &right 
) 
{ 
    float filmAspectRatio = filmApertureWidth / filmApertureHeight; 
    float deviceAspectRatio = imageWidth / (float)imageHeight; 
 
    top = ((filmApertureHeight * inchToMm / 2) / focalLength) * nearClippingPLane; 
    right = ((filmApertureWidth * inchToMm / 2) / focalLength) * nearClippingPLane; 
 
    // field of view (horizontal)
    float fov = 2 * 180 / M_PI * atan((filmApertureWidth * inchToMm / 2) / focalLength); 
    std::cerr << "Field of view " << fov << std::endl; 
 
    float xscale = 1; 
    float yscale = 1; 
 
    switch (fitFilm) { 
        default: 
        case kFill: 
            if (filmAspectRatio > deviceAspectRatio) { 
                xscale = deviceAspectRatio / filmAspectRatio; 
            } 
            else { 
                yscale = filmAspectRatio / deviceAspectRatio; 
            } 
            break; 
        case kOverscan: 
            if (filmAspectRatio > deviceAspectRatio) { 
                yscale = filmAspectRatio / deviceAspectRatio; 
            } 
            else { 
                xscale = deviceAspectRatio / filmAspectRatio; 
            } 
            break; 
    } 
 
    right *= xscale; 
    top *= yscale; 
 
    bottom = -top; 
    left = -right; 
}

void convertToRaster( 
    const glm::vec4 &vertexWorld, 
    const glm::mat4 &worldToCamera, 
    const float &l, 
    const float &r, 
    const float &t, 
    const float &b, 
    const float &near, 
    const uint32_t &imageWidth, 
    const uint32_t &imageHeight, 
    glm::vec3 &vertexRaster 
) 
{ 
    glm::vec3 vertexCamera; 
 	
 	vertexCamera = worldToCamera * vertexWorld;
 
    // convert to screen space
    glm::vec2 vertexScreen; 
    vertexScreen.x = near * vertexCamera.x / -vertexCamera.z; 
    vertexScreen.y = near * vertexCamera.y / -vertexCamera.z; 
 
    // now convert point from screen space to NDC space (in range [-1,1])
    glm::vec2 vertexNDC; 
    vertexNDC.x = 2 * vertexScreen.x / (r - l) - (r + l) / (r - l); 
    vertexNDC.y = 2 * vertexScreen.y / (t - b) - (t + b) / (t - b); 
 
    // convert to raster space
    vertexRaster.x = (vertexNDC.x + 1) / 2 * imageWidth; 
    // in raster space y is down so invert direction
    vertexRaster.y = (1 - vertexNDC.y) / 2 * imageHeight; 
    vertexRaster.z = -vertexCamera.z; 
} 
 
float min3(const float &a, const float &b, const float &c) 
{ return std::min(a, std::min(b, c)); } 
 
float max3(const float &a, const float &b, const float &c) 
{ return std::max(a, std::max(b, c)); } 
 
float edgeFunction(const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c) 
{ return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]); } 
 
const uint32_t imageWidth = 640; 
const uint32_t imageHeight = 480; 
float indexes[16] = {0.707107, -0.331295, 0.624695, 0, 0, 0.883452, 0.468521, 0, -0.707107, -0.331295, 0.624695, 0, -1.63871, -5.747777, -40.400412, 1};
const glm::mat4 worldToCamera = glm::make_mat4(indexes);
 
const uint32_t ntris = 3156; 
const float nearClippingPLane = 1; 
const float farClippingPLane = 1000; 
float focalLength = 20; // in mm 
// 35mm Full Aperture in inches
float filmApertureWidth = 0.980; 
float filmApertureHeight = 0.735;

int main(){
	
	glm::mat4 cameraToWorld = glm::inverse(worldToCamera); 
	 
	// compute screen coordinates
	float t, b, l, r; 
	 
	computeScreenCoordinates( 
	    filmApertureWidth, filmApertureHeight, 
	    imageWidth, imageHeight, 
	    kOverscan, 
	    nearClippingPLane, 
	    focalLength, 
	    t, b, l, r); 

	unsigned char * frameBuffer = new unsigned char[imageHeight*imageWidth*3];
    for (uint32_t i = 0; i < 3*imageWidth * imageHeight; i+=3){
    	frameBuffer[i] = 255; frameBuffer[i+1] = 255; frameBuffer[i+2] = 255; 
    }
	float *depthBuffer = new float[imageWidth * imageHeight]; 
	for (uint32_t i = 0; i < imageWidth * imageHeight; ++i) depthBuffer[i] = farClippingPLane; 
	 
	auto t_start = std::chrono::high_resolution_clock::now();

	for (uint32_t i = 0; i < ntris; ++i) { 

		const glm::vec4 &v0 = glm::vec4(vertices[nvertices[i*3]].x, vertices[nvertices[i*3]].y, vertices[nvertices[i*3]].z, 1);
	    const glm::vec4 &v1 = glm::vec4(vertices[nvertices[i*3+1]].x, vertices[nvertices[i*3+1]].y, vertices[nvertices[i*3+1]].z, 1);
	    const glm::vec4 &v2 = glm::vec4(vertices[nvertices[i*3+2]].x, vertices[nvertices[i*3+2]].y, vertices[nvertices[i*3+2]].z, 1);

	    // Convert vertices triangles to raster space
	    glm::vec3 v0Raster, v1Raster, v2Raster;
	    convertToRaster(v0, worldToCamera, l, r, t, b, nearClippingPLane, imageWidth, imageHeight, v0Raster); 
	    convertToRaster(v1, worldToCamera, l, r, t, b, nearClippingPLane, imageWidth, imageHeight, v1Raster); 
	    convertToRaster(v2, worldToCamera, l, r, t, b, nearClippingPLane, imageWidth, imageHeight, v2Raster); 

	    v0Raster.z = 1 / v0Raster.z, 
	    v1Raster.z = 1 / v1Raster.z, 
	    v2Raster.z = 1 / v2Raster.z; 

	    glm::vec2 st0 = st[stindices[i * 3]]; 
	    glm::vec2 st1 = st[stindices[i * 3 + 1]]; 
	    glm::vec2 st2 = st[stindices[i * 3 + 2]]; 
	 
	    st0 *= v0Raster.z, st1 *= v1Raster.z, st2 *= v2Raster.z;

	    float xmin = min3(v0Raster.x, v1Raster.x, v2Raster.x); 
	    float ymin = min3(v0Raster.y, v1Raster.y, v2Raster.y); 
	    float xmax = max3(v0Raster.x, v1Raster.x, v2Raster.x); 
	    float ymax = max3(v0Raster.y, v1Raster.y, v2Raster.y); 

	    // the triangle is out of screen
	    if (xmin > imageWidth - 1 || xmax < 0 || ymin > imageHeight - 1 || ymax < 0) continue;

	    // be careful xmin/xmax/ymin/ymax can be negative. Don't cast to uint32_t
	    uint32_t x0 = std::max(int32_t(0), (int32_t)(std::floor(xmin))); 
	    uint32_t x1 = std::min(int32_t(imageWidth) - 1, (int32_t)(std::floor(xmax))); 
	    uint32_t y0 = std::max(int32_t(0), (int32_t)(std::floor(ymin))); 
	    uint32_t y1 = std::min(int32_t(imageHeight) - 1, (int32_t)(std::floor(ymax))); 
	 
	    float area = edgeFunction(v0Raster, v1Raster, v2Raster);

	    for (uint32_t y = y0; y <= y1; ++y) { 
	        for (uint32_t x = x0; x <= x1; ++x) { 
	            glm::vec3 pixelSample(x + 0.5, y + 0.5, 0); 
	            float w0 = edgeFunction(v1Raster, v2Raster, pixelSample); 
	            float w1 = edgeFunction(v2Raster, v0Raster, pixelSample); 
	            float w2 = edgeFunction(v0Raster, v1Raster, pixelSample); 
	            if (w0 >= 0 && w1 >= 0 && w2 >= 0) { 
	                w0 /= area; 
	                w1 /= area; 
	                w2 /= area; 
	                float oneOverZ = v0Raster.z * w0 + v1Raster.z * w1 + v2Raster.z * w2; 
	                float z = 1 / oneOverZ; 

	                if (z < depthBuffer[y * imageWidth + x]) { 
	                    depthBuffer[y * imageWidth + x] = z; 
	 
	                    glm::vec2 st = st0 * w0 + st1 * w1 + st2 * w2; 
	 
	                    st *= z; 

	                    glm::vec3 v0Cam, v1Cam, v2Cam; 
	                    v0Cam = worldToCamera*v0;
	                    v1Cam = worldToCamera*v1;
	                    v2Cam = worldToCamera*v2;
	 
	                    float px = (v0Cam.x/-v0Cam.z) * w0 + (v1Cam.x/-v1Cam.z) * w1 + (v2Cam.x/-v2Cam.z) * w2; 
	                    float py = (v0Cam.y/-v0Cam.z) * w0 + (v1Cam.y/-v1Cam.z) * w1 + (v2Cam.y/-v2Cam.z) * w2; 
	 
	                    glm::vec3 pt(px * z, py * z, -z); // pt is in camera space 

	                    glm::vec3 n = glm::normalize(glm::cross(v1Cam - v0Cam, v2Cam - v0Cam));
	                    glm::vec3 viewDirection = -pt;
	                    viewDirection = glm::normalize(viewDirection);
	 
	                    float nDotView =  std::max(0.f, glm::dot(n, viewDirection));

	                    const int M = 10; 
	                    float checker = (fmod(st.x * M, 1.0) > 0.5) ^ (fmod(st.y * M, 1.0) < 0.5); 
	                    float c = 0.3 * (1 - checker) + 0.7 * checker; 
	                    nDotView *= c; 
	                    
	                    frameBuffer[y * imageWidth * 3 + x * 3] = nDotView * 255; 
	                    frameBuffer[y * imageWidth * 3 + x * 3 + 1] = nDotView * 255; 
	                    frameBuffer[y * imageWidth * 3 + x * 3 + 2] = nDotView * 255; 
	                }
	            }
	        }
	    }
	}

    auto t_end = std::chrono::high_resolution_clock::now(); 
	auto passedTime = std::chrono::duration<double, std::milli>(t_end - t_start).count(); 
	std::cerr << "Wall passed time:  " << passedTime << " ms" << std::endl;

	std::ofstream ofs; 
    ofs.open("./cow.ppm"); 
    ofs << "P6\n" << imageWidth << " " << imageHeight << "\n255\n";
    ofs.write((char*)frameBuffer, imageWidth * imageWidth * 3); 
    ofs.close(); 
 
    delete [] frameBuffer; 
    delete [] depthBuffer; 

    return 0;
}
