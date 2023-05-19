from django.http import JsonResponse
import subprocess
from django.http import HttpResponse

def simulation_endpoint(request):
    # Access request data
    param1 = request.GET.get('param1')
    param2 = request.GET.get('param2')

    # Execute C++ simulation or related logic
    # ...

    # Format and return the response
    response = {
        #'result': result,
        'message': 'Simulation executed successfully.'
    }
    return JsonResponse(response)
def run_simulation(request):
    # Run the C++ executable
    try:
        subprocess.run(["C:\\Users\\vande\\Programming\\PICplusplus\\build\\Debug\\PIC++Main.exe"])
        return HttpResponse("C++ program executed successfully.")
    except Exception as e:
        return HttpResponse(f"Error executing C++ program: {str(e)}")
