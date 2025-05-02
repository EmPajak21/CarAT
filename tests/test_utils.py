from utils import get_example_data
import pickle


# Test for get_example_data
def test_get_example_data(tmp_path):
    # Create a mock pickle file
    mock_data = {"duplets": {"A", "B"}, "triplets": {"C", "D"}}
    mock_file = tmp_path / "mock_data.pkl"
    with open(mock_file, "wb") as f:
        pickle.dump(mock_data, f)

    # Call the function
    result = get_example_data(file_path=str(mock_file))

    # Assertions
    assert "nodes" in result
    assert result["nodes"] == {"A", "B", "C", "D"}
