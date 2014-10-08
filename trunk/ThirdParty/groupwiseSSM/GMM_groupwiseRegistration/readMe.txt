The input for group-wise registration are two parameters
1. subj, which is a cell array with each element being a struct. The struct consists of vertices and faces of one training set
2. options for the registration which is explained in GMM_groupwiseRegistration.m

To start, go to the directory "data" and load the sample inputs. The run:
[SSM correspondences] = GMM_groupwiseRegistration(subj, gmmOpt);

The SSM is the model built from the training set and the correspondences are the corresponding points on the training set. You will need it later, when we build multi-object pose+shape model.

