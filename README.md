# Hybrid-predictor
# CNN-based hybrid prediction scheme for image reversible data hiding
# This repo contains the official implementation of our paper:
# [**CNN-based hybrid prediction scheme for image reversible data hiding**]
##### Table of Contents  
1. [Introduction](#Introduction)  
2. [Requriements](#Requriements)
## Introduction
Most existing RDH algorithms rely on a single predictor based on the local textures of the image. However, it is observed that exploiting different predictors based on different local textures might achieve better prediction performance. For example, the traditional rhombus predictor (RP) may be more suitable for smooth regions, while complex texture regions may be more appropriate to exploit the CNN-based predictor (CNNP). Therefore, in this paper, a hybrid prediction scheme is proposed that adaptively chooses different predictors based on the local texture complexity to achieve better performance. Experimental results demonstrate that our algorithm could obtain better performance than the state-of-the-art works.
## Requriements

torch==1.6.0+cu101

numpy== 1.21.5

cv2==3.4.2.16

torchvision==0.7.0+cu101

matplotlib== 3.5.3

tqdm==4.26.0

Pillow==9.1.1
