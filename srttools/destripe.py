import numpy as np


def mask_zeros(image, npix_tol=2):
    mask = np.ones(image.shape, dtype=bool)

    good_hor = 0
    for i in range(image.shape[0]):
        line = image[i, :]
        if len(line[line == 0]) > npix_tol:
            mask[i, :] = False
        else:
            good_hor += 1

    good_ver = 0
    for i in range(image.shape[1]):
        line = image[:, i]
        if len(line[line == 0]) > npix_tol:
            mask[:, i] = False
        else:
            good_ver += 1

    masked_image = image[mask].reshape((good_hor, good_ver))
    return masked_image, mask


def clip_and_smooth(img, clip_sigma=3, smooth_window=10, direction=0):
    """
    Examples
    --------
    >>> img = np.zeros((2,2))
    >>> np.all(clip_and_smooth(img) == img)
    True
    >>> img = np.array([[0, 0], [1, 1]])
    >>> np.all(clip_and_smooth(img, direction=0) == img)
    True
    >>> img = np.array([[0, 1], [0, 1]])
    >>> np.all(clip_and_smooth(img, direction=1) == img)
    True
    >>> img = np.array([[1, 1.], [8., 1]])
    >>> np.allclose(clip_and_smooth(img, clip_sigma=1, smooth_window=0),
    ...             [[1, 1], [3.0310889132455352, 1]])
    True
    """
    from scipy.ndimage import gaussian_filter1d
    rms = np.std(img)
    median = np.median(img)
    bad = img - median > clip_sigma * rms
    img[bad] = clip_sigma * rms
    bad = median - img > clip_sigma * rms
    img[bad] = - clip_sigma * rms
    img = gaussian_filter1d(img, smooth_window, axis=np.logical_not(direction))
    return img


def basket_weaving(img_hor, img_ver, clip_sigma=3, niter_max=4):
    """Basket-Weaving algorithm from Mueller et al. 1707.05573v6."""
    it = 1
    img_hor = np.copy(img_hor)
    img_ver = np.copy(img_ver)
    width = np.max(img_hor.shape)

    while it < niter_max:
        window = width // 2**it
        if window < 4:
            break
        diff = img_hor - img_ver
        diff = clip_and_smooth(diff, clip_sigma=clip_sigma,
                               smooth_window=window, direction=0)

        img_hor = img_hor - diff

        diff1 = img_ver - img_hor
        diff1 = clip_and_smooth(diff1, clip_sigma=clip_sigma,
                                smooth_window=window, direction=1)

        img_ver = img_ver - diff1
        it += 1

    return (img_ver + img_hor) / 2


def destripe_wrapper(image_hor, image_ver, alg='basket-weaving',
                     niter=4):
    image_mean = (image_hor + image_ver) / 2
    masked_image, mask = mask_zeros(image_mean)

    image_mean[mask] = \
        basket_weaving(image_hor[mask].reshape(masked_image.shape),
                       image_ver[mask].reshape(masked_image.shape),
                       niter_max=niter).flatten()
    if alg == 'basket-weaving':
        return image_mean
