---
title: Convolutional neural network for handwritten digit recognition
layout: project
name: Image Recognition of Handwritten Digits
---

[Convolutional neural network](https://towardsdatascience.com/a-comprehensive-guide-to-convolutional-neural-networks-the-eli5-way-3bd2b1164a53), invented by [Yann LeCun](http://yann.lecun.com) *et al*. and later improved by many others, has been effective on image recognition tasks.

In this tutorial I will show how to use [Keras](https://keras.io) API to implement a simple Convolutional Neural Netwrok to do handwritten digit recognition, which could achieve an accuracy as high as 95%.

### We first load the MNIST data downloaded from the internet


```python
import struct
import numpy as np
def read_idx(filename):
    with open(filename, 'rb') as f:
        zero, data_type, dims = struct.unpack('>HBB', f.read(4))
        shape = tuple(struct.unpack('>I', f.read(4))[0] for d in range(dims))
        return np.fromstring(f.read(), dtype=np.uint8).reshape(shape)

train_images_filename = "train-images-idx3-ubyte"
train_labels_filename = "train-labels-idx1-ubyte"
test_images_filename = "t10k-images-idx3-ubyte"
test_labels_filename = "t10k-labels-idx1-ubyte"
train_images = read_idx(train_images_filename)
train_labels = read_idx(train_labels_filename)
test_images = read_idx(test_images_filename)
test_labels = read_idx(test_labels_filename)
```

    /usr/local/lib/python3.7/site-packages/ipykernel_launcher.py:7: DeprecationWarning: The binary mode of fromstring is deprecated, as it behaves surprisingly on unicode inputs. Use frombuffer instead
      import sys


### Now we reshape the data to make it appropriate for ConvNet.


```python
from keras.utils import np_utils
print(train_images.shape, train_labels.shape, test_images.shape, test_labels.shape)
train_images = np.expand_dims(train_images, axis = -1)
test_images = np.expand_dims(test_images, axis = -1)
train_labels = np_utils.to_categorical(train_labels)
print(train_images.shape, train_labels.shape, test_images.shape, test_labels.shape)
```

    Using TensorFlow backend.


    (60000, 28, 28) (60000,) (10000, 28, 28) (10000,)
    (60000, 28, 28, 1) (60000, 10) (10000, 28, 28, 1) (10000,)


### Before we proceed, we plot one example train data out.


```python
import matplotlib.pyplot as plt
%matplotlib inline
plt.imshow(train_images[0,:,:,0])
print(train_labels[0].argmax())
```

    5



![png](/assets/img/sampleImg.png)


### Now we build a sequential ConvNet using Keras API.
In our sequential mode, we have 3 Convolutional layers, 2 maxpooling layers, followed by 2 dense layers and a final softmax activation. Since the problem is multi-cateogorical imge recogonition, we used **cateogorical cross entropy** as our loss function. We used **adam** optimizer for parameter optimization.
1 **Epoch** of optimization is good enough to acheive a 95% accuracy.


```python
from keras.layers import Conv2D, Dense, MaxPooling2D, Flatten, BatchNormalization
from keras.models import Sequential
from keras.utils import np_utils


mySimpleModel = Sequential()
mySimpleModel.add(Conv2D(32, kernel_size = (3,3), activation='relu', padding = 'same', input_shape=(28,28,1)))
mySimpleModel.add(MaxPooling2D(pool_size=(2,2)))
mySimpleModel.add(BatchNormalization())
mySimpleModel.add(Conv2D(64, kernel_size = (3,3), activation='relu', padding = 'same'))
mySimpleModel.add(MaxPooling2D(pool_size=(2,2)))
mySimpleModel.add(BatchNormalization())
mySimpleModel.add(Conv2D(128, kernel_size = (3,3), activation='relu', padding = 'same'))
mySimpleModel.add(Flatten())
mySimpleModel.add(BatchNormalization())
mySimpleModel.add(Dense(64))
mySimpleModel.add(Dense(10, activation = 'softmax'))

mySimpleModel.compile(loss='categorical_crossentropy', optimizer = 'adam', metrics= ['accuracy'])

mySimpleModel.fit(train_images, train_labels, epochs = 1, verbose = 1)
```

    WARNING:tensorflow:From /usr/local/lib/python3.7/site-packages/tensorflow/python/framework/op_def_library.py:263: colocate_with (from tensorflow.python.framework.ops) is deprecated and will be removed in a future version.
    Instructions for updating:
    Colocations handled automatically by placer.
    WARNING:tensorflow:From /usr/local/lib/python3.7/site-packages/tensorflow/python/ops/math_ops.py:3066: to_int32 (from tensorflow.python.ops.math_ops) is deprecated and will be removed in a future version.
    Instructions for updating:
    Use tf.cast instead.
    Epoch 1/1
    60000/60000 [==============================] - 156s 3ms/step - loss: 0.4089 - acc: 0.9536





    <keras.callbacks.History at 0x13242fc18>



### Now we use the trained model to predict some test images and calculate the ratio of misclassified images and correctly classified images.


```python
test_pred = mySimpleModel.predict(test_images)
test_pred_cat = np.argmax(test_pred, axis = 1)
mismatch_labels = test_pred_cat[test_pred_cat != test_labels]
mismatch_images = test_images[test_pred_cat != test_labels]
mismatch_test = test_labels[test_pred_cat != test_labels]
mismatch = zip(mismatch_labels, mismatch_images, mismatch_test)
match_labels = test_pred_cat[test_pred_cat == test_labels]
match_images = test_images[test_pred_cat == test_labels]
match_test = test_labels[test_pred_cat == test_labels]
match = zip(match_labels, match_images, match_test)
print(len(mismatch_test)*1.0/len(test_labels))
print(len(match_test)/len(test_labels))
```

    0.0292
    0.9708


### Now we randomly choose 10 images from the misclassified set and plot them.


```python
#random choose 10 from mismatches to display
choices = np.random.choice(range(len(mismatch_labels)),size=10, replace=False)
fig, axs = plt.subplots(5,2,sharex=True,figsize=(20,40))
for i in range(len(choices)):
    index = choices[i]
    ax = axs[i//2][i%2]
    img = mismatch_images[index]
    mislabel = mismatch_labels[index]
    realLabel = mismatch_test[index]
    ax.imshow(img[:,:,0])
    ax.set_title('mislabel:'+str(mislabel)+','+'realLabel:'+str(realLabel))
plt.show()
```


![png](/assets/img/misclassified.png)


### Now we randomly choose 10 images from the correctly classified set and plot them.


```python
#random choose 10 from mismatches to display
choices = np.random.choice(range(len(match_labels)),size=10, replace=False)
fig, axs = plt.subplots(5,2,sharex=True,figsize=(20,40))
for i in range(len(choices)):
    index = choices[i]
    ax = axs[i//2][i%2]
    img = match_images[index]
    label = match_labels[index]
    realLabel = match_test[index]
    ax.imshow(img[:,:,0])
    ax.set_title('mislabel:'+str(label)+','+'realLabel:'+str(realLabel))
plt.show()
```


![png](/assets/img/correctclassified.png)