
import glob
import torch
import numpy as np
import os
import cv2
import matlab
import matlab.engine
import utils
import math
from model.predict_model import PredictModel


def PsnrC(img1, img2):
    # calculate the PSNR
    img1 = np.array(img1, dtype=np.float)
    img2 = np.array(img2, dtype=np.float)
    mse = np.mean((img1 - img2) ** 2)
    if mse < 1.0e-10:
        return 100
    return 10 * math.log10(255.0 ** 2 / mse)


def main():
    engine = matlab.engine.start_matlab()  # start Matlab process
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--model-pth', '-model', default=r'.\model_parameter\model_state.pth', type=str,
                        help='The place of the model parameter.')  # model
    parser.add_argument('--image', type=str, default=None)  # The image to be embedded（single）
    parser.add_argument('--images_dir', type=str,
                        default=r'D:\pengyi\master_student\work_CNN\xwcnnmhs_final\test_images_dir')  # The image to be embedded

    parser.add_argument('--save_dir', type=str,
                        default=r'D:\pengyi\master_student\work_CNN\xwcnnmhs_final\stego_images_dir')  # stego image

    parser.add_argument('--secret', type=int, default=0.01)  # secret

    args = parser.parse_args()
    if args.image is not None:
        files_list = [args.image]
    elif args.images_dir is not None:
        files_list = glob.glob(args.images_dir + '/*')
    else:
        print('Missing input image')
        return

    # model_path = args.model
    m = args.secret

    width = 512
    height = 512

    size = (width, height)

    img_psnr1 = []

    device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
    model = PredictModel(device)

    utils.load_model(args.model_pth, model)

    for i in range(0, len(files_list)):
        print(i)

        img_path = args.images_dir
        filename = os.path.join(img_path, files_list[i])
        img = cv2.imread(filename)
        img_resize = cv2.resize(img, size, interpolation=cv2.INTER_CUBIC)
        img_gray = cv2.cvtColor(img_resize, cv2.COLOR_BGR2GRAY)
        img_gray = np.array(img_gray, dtype=np.float)

        stego_img = np.zeros(img_gray.shape)
        stego_img = utils.cnn_mhs(img_gray, m, device, model, engine)

        img_psnr1.append(np.round(PsnrC(stego_img, img_gray), 5))

    print(img_psnr1)
    print('CNNP psnr1 = {}'.format(np.round(np.mean(img_psnr1), 5)))

    engine.exit()


if __name__ == "__main__":
    main()
