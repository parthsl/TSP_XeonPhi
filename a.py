import numpy as np
import cv2

cap = cv2.VideoCapture(0)
fp = None
while(True):
    # Capture frame-by-frame
    ret, frame = cap.read()
    # Our operations on the frame come here
    #gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    gray = frame
    kernel = np.ones((3,3), np.float32)/9
    #cv2.filter2D(frame, -1, kernel, dst=gray)
    # Display the resulting frame
    cv2.imshow('frame',gray)
    waitk = cv2.waitKey(1)
    if (waitk == ord('q')):
        break
    if waitk == ord('c'):
        fp = gray
        print(fp,"dsfdfdf")
# When everything done, release the capture

while(True):
    ret, frame = cap.read()
    # Our operations on the frame come here
    #gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    gray = frame
    kernel = np.ones((3,3), np.float32)/9
    cv2.filter2D(frame, -1, fp, dst=gray)
    # Display the resulting frame
    cv2.imshow('frame',gray)
    waitk = cv2.waitKey(1)
    if (waitk == ord('z')):
        break


cap.release()
cv2.destroyAllWindows()
