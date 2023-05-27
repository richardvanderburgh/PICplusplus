from django.http import JsonResponse
import subprocess
from django.http import HttpResponse
import json

def run_simulation(request):
    if request.method == 'GET':
        param1 = request.GET.get('param1')
        param2 = request.GET.get('param2')
    else:
        param1 = request.POST.get('param1')
        param2 = request.POST.get('param2')

    # Run the C++ executable
    try:
        executablePath = "C:\\Users\\vande\\Programming\\PICplusplus\\build\\bin\\PIC++Main.exe"
        result = subprocess.run([executablePath, param1, param2], capture_output=True, text=True)
        output = result.stdout.strip()

        return JsonResponse(output, safe=False)
        #return HttpResponse("C++ program executed successfully.")
    except Exception as e:
        return HttpResponse(f"Error executing C++ program: {str(e)}")
