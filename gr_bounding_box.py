__author__ = 'christopher'

import numpy as np

def gr_intensity_error(r, q, dq, dsq):
    error = np.zeros(r.shape)
    for i in range(len(q)):
        error[:] += (q[i]*np.sin(q[i]*r[:])*dq[i])**2*dsq[i]
    error *= 4/np.pi**2
    return error

def get_total_bounding_box(qmax, qmin, rmax, rmin, dq, dsq):
    q = np.arange(qmin, qmax, dq)
    dq = np.ones(q.shape)*dq
    r = np.arange(rmin, rmax, np.pi/qmax)
    return np.pi**3/qmax**2*np.sum(gr_intensity_error(r, q, dq, dsq))


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from Mask.IO import load_chi_file
    file_loc = '/mnt/bulk-data/research_data/Brown_beamtime/Beamtime/Mar2014/7_7_7_NiPd_Ambient/NiPd_LongRun/7_7_7_NiPd_AsIs_300K-00020_sum_qspace.chi'
    data, filenname = load_chi_file(file_loc, 29)
    print data.shape
    # plt.plot(data[:,0], data[:,1], label='data')
    # plt.plot(data[:,0], data[:,2], label='var')
    # plt.legend()
    # plt.show()
    dsq = data[:, -1]/12.**4
    rmax = 40.
    rmin = 0.
    dq = .1
    qmin = 0.
    # dsq = np.ones(25/.1)
    errors = []
    q = np.arange(0, 25, .1)
    ierror = []
    rerror = []
    for qmax in q:
        errors.append(get_total_bounding_box(qmax, qmin, rmax, rmin, dq, dsq))
        ierror.append(errors[-1]*qmax**2/np.pi**2)
        rerror.append(np.pi**2/qmax**2)


    plt.plot(q, errors, label='total')
    plt.plot(q, ierror, label='intensity')
    plt.plot(q, rerror, label='r')
    plt.legend()
    # plt.plot([q[0], q[-1]], [errors[0], errors[-1]])
    plt.show()
