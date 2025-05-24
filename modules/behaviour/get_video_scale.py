import cv2
import matplotlib.pyplot as plt

def get_scale_from_frame(video_path, real_world_distance_cm, frame_number=10):
    """
    Lets the user click two points in a video frame to calculate the cm/px scale.
    
    Parameters:
        video_path (str): Path to the video file.
        real_world_distance_cm (float): Known distance between two points in cm.
        frame_number (int): Frame to grab from the video (default: 0).
        
    Returns:
        scale (float): Conversion factor in cm per pixel.
    """
    # Load the video
    cap = cv2.VideoCapture(video_path)
    cap.set(cv2.CAP_PROP_POS_FRAMES, frame_number)
    success, frame = cap.read()
    cap.release()
    
    if not success:
        raise ValueError(f"Failed to read frame {frame_number} from {video_path}")
    
    # Convert BGR (OpenCV) to RGB (Matplotlib)
    frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    
    # Show the frame and get two points
    plt.imshow(frame_rgb)
    plt.title("Click two points with known real-world distance")
    points = plt.ginput(2, timeout=0)
    plt.close()
    
    if len(points) != 2:
        raise ValueError("You must click exactly two points.")
    
    # Compute pixel distance
    p1, p2 = points
    pixel_distance = ((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)**0.5
    
    # Calculate scale
    scale = real_world_distance_cm / pixel_distance
    print(f"Pixel distance: {pixel_distance:.2f} px")
    print(f"Scale: {scale:.4f} cm/px")
    
    return scale

if __name__ == "__main__":
    video_folder = r'E:\FiberPhotometry\202504_OptoFluidACh\DLC_Projects\FiberMEC_EPM-Alice-2025-05-14\videos_original'
    video_name = '765_0_reduced'
    video_path = f'{video_folder}\{video_name}.avi'
    output_json = f'{video_name}_epm_boundaries.json'
    real_world_distance_cm = 35

    scale = get_scale_from_frame(video_path, real_world_distance_cm, frame_number=10)
    print(f"Scale in {video_name} : {scale}px/cm")